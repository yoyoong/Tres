import pandas as pd
from rdkit import Chem
from rdkit import RDLogger


def toStandardSmile(train_df):
    RDLogger.DisableLog('rdApp.*')  # 屏蔽RDKit的warning
    for index, row in train_df.iterrows():
        try:
            mol = Chem.MolFromSmiles(row['SMILES'])
            new_smiles = Chem.MolToSmiles(mol)
            train_df.loc[index, 'SMILES'] = new_smiles
        except:
            # 若转化失败，则认为原始smile不合法，删除该数据
            train_df.drop(index, inplace=True)
    return train_df

def removeDuplicate(train_df):
    duplicate_rows = train_df[train_df.duplicated('SMILES', keep=False)]
    for smiles, group in duplicate_rows.groupby('SMILES'):
        if len(group.drop_duplicates(subset=['label'])) == 1:
            train_df.drop(index=group.index[1:], inplace=True)
        else:
            train_df.drop(index=group.index, inplace=True)
    return train_df

if __name__ == '__main__':
    # # 导入训练集
    # train_df = pd.read_csv('data/train.csv')
    # print(f'len of train_df is {len(train_df)}')
    #
    # # 将smiles转化为rdkit标准smiles
    # train_df = toStandardSmile(train_df)
    #
    # # 去除重复值
    # train_df = removeDuplicate(train_df)
    # print(f'len of train_df after pre-pocess is {len(train_df)}')
    #
    # # 保存处理好的csv文件
    # train_df.to_csv('train_preprocessed.csv', index=0)

    # 读取预处理好的训练数据
    train_df = pd.read_csv('train_preprocessed.csv')

    # 读取测试集
    test_df = pd.read_csv('data/test_nolabel.csv')
    smiles_list = test_df["SMILES"].tolist()

    torch.load('./model/pretrained_synpg.pt')





