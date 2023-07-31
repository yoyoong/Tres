import pandas as pd
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import AllChem
import seaborn as sns
import matplotlib.pyplot as plt
import pickle as pkl
import pandas as pd
from threading import Thread, Lock
from pahelix.utils.compound_tools import mol_to_geognn_graph_data_MMFF3d
import paddle as pdl
import paddle.nn as nn
from pahelix.model_zoo.gem_model import GeoGNNModel
import json
import pgl
from pahelix.datasets.inmemory_dataset import InMemoryDataset
from sklearn.model_selection import train_test_split
from paddle import optimizer
import numpy as np
import paddle.nn.functional as F
from sklearn.metrics import average_precision_score, roc_auc_score, accuracy_score, precision_score, recall_score, \
    f1_score, confusion_matrix
import warnings, os

warnings.filterwarnings('ignore')
os.environ['FLAGS_cudnn_deterministic'] = 'True'

pdl.seed(123)
np.random.seed(2021)

# 导入训练集
train_df = pd.read_csv('data/train.csv')
print(f'len of train_df is {len(train_df)}')

# 将smiles转化为rdkit标准smiles
RDLogger.DisableLog('rdApp.*')  # 屏蔽RDKit的warning
for index, row in train_df.iterrows():
    try:
        mol = Chem.MolFromSmiles(row['SMILES'])
        new_smiles = Chem.MolToSmiles(mol)
        train_df.loc[index, 'SMILES'] = new_smiles
    except:
        # 若转化失败，则认为原始smile不合法，删除该数据
        train_df.drop(index, inplace=True)

# 去除重复值
duplicate_rows = train_df[train_df.duplicated('SMILES', keep=False)]
for smiles, group in duplicate_rows.groupby('SMILES'):
    if len(group.drop_duplicates(subset=['label'])) == 1:
        train_df.drop(index=group.index[1:], inplace=True)
    else:
        train_df.drop(index=group.index, inplace=True)
print(f'len of train_df after pre-pocess is {len(train_df)}')

# 观察label的分布
# sns.set()
# sns.histplot(train_df['label'])
# plt.show()

# 保存处理好的csv文件
train_df.to_csv('train_preprocessed.csv', index=0)

# 将smiles列表保存为smiles_list.pkl文件
train_df = pd.read_csv('train_preprocessed.csv')  # 读取预处理好的训练数据
smiles_list = train_df["SMILES"].tolist()
pkl.dump(smiles_list, open('GEM/work/train_smiles_list.pkl', 'wb'))

# 测试集
test_df = pd.read_csv('data/test_nolabel.csv')
smiles_list = test_df["SMILES"].tolist()
pkl.dump(smiles_list, open('GEM/work/test_smiles_list.pkl', 'wb'))

# 使用分子力场将smiles转化为3d分子图，并保存为smiles_to_graph_dict.pkl文件
mutex = Lock()  # 互斥锁，防止多个线程同时修改某一文件或某一全局变量，引发未知错误

def calculate_3D_structure_(smiles_list):
    n = len(smiles_list)
    global p
    index = 0
    while True:
        mutex.acquire()  # 获取锁
        if p >= n:
            mutex.release()
            break
        index = p  # p指针指向的位置为当前线程要处理的smiles
        smiles = smiles_list[index]
        print(index, ':', round(index / n * 100, 2), '%', smiles)
        p += 1  # 修改全局变量p
        mutex.release()  # 释放锁
        try:
            molecule = AllChem.MolFromSmiles(smiles)
            molecule_graph = mol_to_geognn_graph_data_MMFF3d(molecule)  # 根据分子力场生成3d分子图
        except:
            print("Invalid smiles!", smiles)
            mutex.acquire()
            with open('GEM/work/invalid_smiles.txt', 'a') as f:
                # 生成失败的smiles写入txt文件保存在该目录下
                f.write(str(smiles) + '\n')
            mutex.release()
            continue

        global smiles_to_graph_dict
        mutex.acquire()  # 获取锁
        smiles_to_graph_dict[smiles] = molecule_graph
        mutex.release()  # 释放锁

# for mode in ['train', 'test']:
#     if mode == 'train':
#         smiles_list = train_df["SMILES"].tolist()
#     else:
#         smiles_list = test_df["SMILES"].tolist()
#     global smiles_to_graph_dict
#     smiles_to_graph_dict = {}
#     global p              # p为全局指针，指向即将要处理的smiles
#     p = 0
#     thread_count = 12      # 线程数。一般根据当前运行环境下cpu的核数来设定
#     threads = []
#     for i in range(thread_count):
#         threads.append(Thread(target=calculate_3D_structure_, args=(smiles_list, )))
#     for t in threads:
#         t.start()
#     for t in threads:
#         t.join()
#     pkl.dump(smiles_to_graph_dict, open(f'work/{mode}_smiles_to_graph_dict.pkl', 'wb'))
#     print(f'{mode} is Done!')


# 将smiles、graph、label构建成一个列表，并保存为data_list.pkl文件，该文件为GEM读取的数据文件
train_smiles_to_graph_dict = pkl.load(open(f'GEM/work/train_smiles_to_graph_dict.pkl', 'rb'))
test_smiles_to_graph_dict = pkl.load(open(f'GEM/work/test_smiles_to_graph_dict.pkl', 'rb'))
train_data_list = []
test_data_list = []

for index, row in train_df.iterrows():
    smiles = row["SMILES"]
    label = row["label"]
    if smiles not in train_smiles_to_graph_dict:
        continue
    data_item = {
        "smiles": smiles,
        "graph": train_smiles_to_graph_dict[smiles],
        "label": label,
    }
    train_data_list.append(data_item)

for index, row in test_df.iterrows():
    smiles = row["SMILES"]
    if smiles not in test_smiles_to_graph_dict:
        continue
    data_item = {
        "smiles": smiles,
        "graph": test_smiles_to_graph_dict[smiles],
        'label': 0
    }
    test_data_list.append(data_item)

pkl.dump(train_data_list, open('GEM/work/train_data_list.pkl', 'wb'))
pkl.dump(test_data_list, open('GEM/work/test_data_list.pkl', 'wb'))


class ADMET(nn.Layer):
    def __init__(self):
        super(ADMET, self).__init__()
        compound_encoder_config = json.load(open('GEM/model_configs/geognn_l8.json', 'r'))
        self.encoder = GeoGNNModel(compound_encoder_config)
        self.encoder.set_state_dict(pdl.load("GEM/weight/class.pdparams"))
        # GEM编码器输出的图特征为32维向量, 因此mlp的输入维度为32
        self.mlp = nn.Sequential(
            nn.Linear(32, 32, weight_attr=nn.initializer.KaimingNormal()),
            nn.ReLU(),
            nn.Linear(32, 32, weight_attr=nn.initializer.KaimingNormal()),
            nn.ReLU(),
            nn.Linear(32, 32, weight_attr=nn.initializer.KaimingNormal()),
            nn.ReLU(),
            nn.Linear(32, 2, weight_attr=nn.initializer.KaimingNormal()),
        )

    def forward(self, atom_bond_graph, bond_angle_graph):
        node_repr, edge_repr, graph_repr = self.encoder(atom_bond_graph.tensor(), bond_angle_graph.tensor())
        return self.mlp(graph_repr)


def collate_fn(data_batch):
    """
    Dataloader中的数据处理函数：该函数输入一个batch的数据, 返回一个batch的(atom_bond_graph, bond_angle_graph, label)
    """
    atom_names = ["atomic_num", "formal_charge", "degree", "chiral_tag", "total_numHs", "is_aromatic", "hybridization"]
    bond_names = ["bond_dir", "bond_type", "is_in_ring"]
    bond_float_names = ["bond_length"]
    bond_angle_float_names = ["bond_angle"]

    atom_bond_graph_list = []  # 原子-键特征图
    bond_angle_graph_list = []  # 键-键角特征图
    label_list = []  # label

    for data_item in data_batch:
        graph = data_item['graph']
        ab_g = pgl.Graph(
            num_nodes=len(graph[atom_names[0]]),
            edges=graph['edges'],
            node_feat={name: graph[name] for name in atom_names},
            edge_feat={name: graph[name] for name in bond_names + bond_float_names})
        ba_g = pgl.Graph(
            num_nodes=len(graph['edges']),
            edges=graph['BondAngleGraph_edges'],
            node_feat={name: graph[name] for name in bond_names + bond_float_names},
            edge_feat={name: graph[name] for name in bond_angle_float_names})
        atom_bond_graph_list.append(ab_g)
        bond_angle_graph_list.append(ba_g)
        label_list.append(data_item['label'])

    atom_bond_graph = pgl.Graph.batch(atom_bond_graph_list)
    bond_angle_graph = pgl.Graph.batch(bond_angle_graph_list)

    return atom_bond_graph, bond_angle_graph, np.array(label_list, dtype=np.float32)


def get_data_loader(mode, batch_size):
    if mode == 'train':
        # 训练模式下将train_data_list划分训练集和验证集，返回对应的DataLoader
        data_list = pkl.load(open(f'GEM/work/train_data_list.pkl', 'rb'))  # 读取data_list

        train_data_list, valid_data_list = train_test_split(data_list, test_size=0.2, random_state=42)
        print(f'train: {len(train_data_list)}, valid: {len(valid_data_list)}')

        train_dataset = InMemoryDataset(train_data_list)
        valid_dataset = InMemoryDataset(valid_data_list)
        train_data_loader = train_dataset.get_data_loader(
            batch_size=batch_size, num_workers=1, shuffle=True, collate_fn=collate_fn)
        valid_data_loader = valid_dataset.get_data_loader(
            batch_size=batch_size, num_workers=1, shuffle=True, collate_fn=collate_fn)
        return train_data_loader, valid_data_loader

    elif mode == 'test':
        # 推理模式下直接读取test_data_list, 返回test_data_loader
        data_list = pkl.load(open(f'GEM/work/test_data_list.pkl', 'rb'))

        test_dataset = InMemoryDataset(data_list)
        test_data_loader = test_dataset.get_data_loader(
            batch_size=batch_size, num_workers=1, shuffle=False, collate_fn=collate_fn)
        return test_data_loader


def trial(model_version, model, batch_size, criterion, scheduler, opt):
    # 创建dataloader
    train_data_loader, valid_data_loader = get_data_loader(mode='train', batch_size=batch_size)

    current_best_metric = -1e10
    max_bearable_epoch = 50  # 设置早停的轮数为50，若连续50轮内验证集的评价指标没有提升，则停止训练
    current_best_epoch = 0

    for epoch in range(1000):  # 设置最多训练1000轮
        train_label_true = pdl.to_tensor([], dtype=pdl.float32, place=pdl.CUDAPlace(0))
        train_label_predict = pdl.to_tensor([], dtype=pdl.float32, place=pdl.CUDAPlace(0))

        model.train()
        for (atom_bond_graph, bond_angle_graph, label_true_batch) in train_data_loader:
            label_predict_batch = model(atom_bond_graph, bond_angle_graph)
            label_true_batch = pdl.to_tensor(label_true_batch, dtype=pdl.int64, place=pdl.CUDAPlace(0))
            train_label_true = pdl.concat((train_label_true, label_true_batch.detach()), axis=0)
            label_predict_batch_softmax = F.softmax(label_predict_batch)
            train_label_predict = pdl.concat((train_label_predict, label_predict_batch_softmax.detach()), axis=0)

            loss = criterion(label_predict_batch, label_true_batch)
            loss.backward()  # 反向传播
            opt.step()  # 更新参数
            opt.clear_grad()
        scheduler.step()  # 更新学习率

        valid_label_true = pdl.to_tensor([], dtype=pdl.float32, place=pdl.CUDAPlace(0))
        valid_label_predict = pdl.to_tensor([], dtype=pdl.float32, place=pdl.CUDAPlace(0))
        model.eval()
        for (atom_bond_graph, bond_angle_graph, label_true_batch) in valid_data_loader:
            label_predict_batch = model(atom_bond_graph, bond_angle_graph)
            label_true_batch = pdl.to_tensor(label_true_batch, dtype=pdl.int64, place=pdl.CUDAPlace(0))
            valid_label_true = pdl.concat((valid_label_true, label_true_batch.detach()), axis=0)
            label_predict_batch_softmax = F.softmax(label_predict_batch)
            valid_label_predict = pdl.concat((valid_label_predict, label_predict_batch_softmax.detach()), axis=0)

        # 评估模型在训练集、验证集的表现
        metric_train = evaluate(train_label_true, train_label_predict)
        metric_valid = evaluate(valid_label_true, valid_label_predict)

        score = round((metric_valid['ap'] + metric_valid['auc']) / 2, 4)

        if score > current_best_metric:
            # 保存score最大时的模型权重
            current_best_metric = score
            current_best_epoch = epoch
            pdl.save(model.state_dict(), "trainedWeight/" + model_version + ".pkl")
        print('current_epoch', epoch, 'current_best_epoch', current_best_epoch, 'current_best_metric',
              current_best_metric)
        if epoch > current_best_epoch + max_bearable_epoch:
            break

    return;


def evaluate(label_true, label_predict):
    """评估模型"""
    y_pred = label_predict[:, 1].cpu().numpy()
    y_true = label_true.cpu().numpy()

    ap = round(average_precision_score(y_true, y_pred), 4)
    auc = round(roc_auc_score(y_true, y_pred), 4)

    y_pred = np.where(y_pred >= 0.5, 1, 0)
    accuracy = round(accuracy_score(y_true, y_pred), 4)
    precision = round(precision_score(y_true, y_pred), 4)
    recall = round(recall_score(y_true, y_pred), 4)
    f1 = round(f1_score(y_true, y_pred), 4)
    confusion_mat = confusion_matrix(y_true, y_pred)

    metric = {'ap': ap, 'auc': auc, 'accuracy': accuracy, 'precision': precision, 'recall': recall, 'f1': f1,
              'confusion_mat': confusion_mat}
    return metric


# 固定随机种子
# SEED = 1024
# pdl.seed(SEED)
# np.random.seed(SEED)
# random.seed(SEED)

model = ADMET()
batch_size = 1024  # batch size
criterion = nn.CrossEntropyLoss()  # 损失函数
scheduler = optimizer.lr.CosineAnnealingDecay(learning_rate=1e-3, T_max=15)  # 余弦退火学习率
opt = optimizer.Adam(scheduler, parameters=model.parameters(), weight_decay=1e-5)  # 优化器
trial(model_version='1', model=model, batch_size=batch_size, criterion=criterion, scheduler=scheduler, opt=opt)


# 将测试集的预测结果保存为result.csv
def test(model_version):
    test_data_loader = get_data_loader(mode='test', batch_size=1024)

    model = ADMET()
    model.set_state_dict(pdl.load("GEM/trainedWeight/" + model_version + ".pkl"))  # 导入训练好的的模型权重
    model.eval()

    all_result = []
    for (atom_bond_graph, bond_angle_graph, label_true_batch) in test_data_loader:
        label_predict_batch = model(atom_bond_graph, bond_angle_graph)
        label_predict_batch = F.softmax(label_predict_batch)
        result = label_predict_batch[:, 1].cpu().numpy().reshape(-1).tolist()
        all_result.extend(result)

    df = pd.read_csv('data/test_nolabel.csv')
    df[f'pred'] = all_result
    df.to_csv('result.csv', index=False)


test(model_version='1')
