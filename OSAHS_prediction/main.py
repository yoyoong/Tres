import numpy as np
import torch
from torch import nn
from torch import optim
from torch.nn import functional as nnFun
import time

def getData():
    train_x = np.load('./A/train_x.npz')['data']  # (21547, 4, 15000)
    train_y = np.load('./A/train_y.npz')['data']  # （21547,)
    test_X = np.load('./A/test_x.npz')['data']  # (1000, 4, 15000)
    return train_x, train_y, test_X

# [m, 4, n] => [m, n, 2, 2]
def vec2matrix(data):
    return data.swapaxes(1, 2).reshape(data.shape[0], data.shape[2], 2, 2)

def try_gpu(i=0):
    """Return gpu(i) if exists, otherwise return cpu()."""
    if torch.cuda.device_count() >= i + 1:
        return torch.device(f'cuda:{i}')
    return torch.device('cpu')

class DatasetLoader(torch.utils.data.Dataset):
    def __init__(self, data, label):
        self.data = data
        self.label = label

    def __getitem__(self, index):
        return self.data[index], self.label[index]

    def __len__(self):
        return len(self.data)

class CNNLSTM(nn.Module):
    def __init__(self, num_classes=4, num_channels=4):
        super(CNNLSTM, self).__init__()

        # 定义CNN部分
        self.cnn = nn.Sequential(
            nn.Conv2d(1, 50, kernel_size=(num_channels, 50)), nn.BatchNorm2d(50),
            nn.MaxPool2d(kernel_size=(1, 50), stride=(1, 25)),
            nn.Conv2d(50, 100, kernel_size=(1, 10)),  nn.BatchNorm2d(100),
            nn.MaxPool2d(kernel_size=(1, 10), stride=(1, 5))
        )

        # 定义LSTM部分
        self.lstm1 = nn.LSTM(11600, 1024, batch_first=True)
        self.lstm2 = nn.LSTM(1024, 128, batch_first=True)

        # 定义全连接层
        self.fc = nn.Linear(128, num_classes)

    def forward(self, x):
        # batch_size, timesteps, C, H, W = x.size()
        batch_size, H, W = x.size()
        C, timesteps = 1, 1

        # CNN部分
        x = x.view(batch_size * timesteps, C, H, W).float()
        x = self.cnn(x)

        # 将特征图转换为LSTM的输入形状
        x = x.view(batch_size, timesteps, -1)

        # LSTM部分
        h0 = torch.zeros(1, batch_size, 1024).to(x.device)
        c0 = torch.zeros(1, batch_size, 1024).to(x.device)
        h1 = torch.zeros(1, batch_size, 128).to(x.device)
        c1 = torch.zeros(1, batch_size, 128).to(x.device)
        lstm_out, _ = self.lstm1(x, (h0, c0))
        lstm_out, _ = self.lstm2(lstm_out, (h1, c1))
        lstm_out = lstm_out[:, -1, :] # 取最后一个时间步的输出

        # 全连接层
        output = self.fc(lstm_out)
        return output

if __name__ == '__main__':
    train_x, train_y, submit_X = getData()
    # trainDataAnalysis(train_x, train_y)
    # train_x = vec2matrix(train_x)

    # params
    batch_size = 100
    num_channels = 4
    num_kernels = 1
    num_timesteps = 15000
    num_epochs = 100
    device = try_gpu()
    model = CNNLSTM().to(device)
    datasetLoader = DatasetLoader(train_x, train_y)
    train_data, test_data = torch.utils.data.random_split(datasetLoader, [18547, 3000])
    train_data = torch.utils.data.DataLoader(train_data, batch_size=batch_size, shuffle=True)
    test_data = torch.utils.data.DataLoader(test_data, batch_size=batch_size, shuffle=False)

    # train
    optimizer = optim.SGD(model.parameters(), lr=0.01, momentum=0.9)
    criterion = nn.CrossEntropyLoss()
    print("Train start!")
    for epoch in range(num_epochs):
        epoch_loss = 0.0
        for batch_index, (x, y) in enumerate(train_data):
            x, y = x.to(device), y.to(device)
            predict = model(x)
            loss = criterion(predict, y)
            epoch_loss += loss.item()
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
        now = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))
        print(now + " Train epoch:" + str(epoch) + "  epoch_mean_loss: " + str(epoch_loss / (batch_index + 1)))
    torch.save(model.state_dict(), 'model1.1.pth')
    print("Train end!")

    # validation
    correct_num = 0
    total_num = 0
    for x, y in test_data:
        x, y = x.to(device), y.to(device)
        predict = model(x)
        predict_value = predict.argmax(axis=1)
        total_num += len(y)
        correct = predict_value.eq(y).sum().float().item()
        correct_num += correct
    print("Predict accuracy:" + str(correct_num / total_num))

