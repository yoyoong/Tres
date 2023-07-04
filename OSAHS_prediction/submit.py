import numpy as np
from main import CNNLSTM, try_gpu, getData
import torch
import csv

device = try_gpu()
model = CNNLSTM()
model.load_state_dict(torch.load('model1.1.pth', map_location='cuda'))
model.to(device)

train_x, train_y, submit_X = getData()

result = []
for x in submit_X:
    x = x.reshape(1, x.shape[0], x.shape[1])
    x = torch.from_numpy(x).to(device)
    predict = model(x)
    predict_value = predict.argmax(axis=1)
    result.append(str(torch.tensor(predict_value, device='cpu').clone().detach().numpy()[0]))

print(result)
with open('submission.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(result)


