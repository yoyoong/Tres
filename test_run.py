import pandas as pd
from lifelines import CoxPHFitter

# 构建示例数据，假设有患者的生存时间、事件发生与一些协变量
data = pd.DataFrame({
    '生存时间': [60, 80, 120, 40, 200, 150, 90, 110],
    '事件发生': [1, 1, 1, 0, 1, 1, 1, 1],
    '年龄': [55, 65, 70, 45, 75, 60, 50, 68],
    '性别': [1, 0, 1, 0, 1, 0, 1, 0],  # 1表示男性，0表示女性
    '癌症类型': ['A', 'B', 'A', 'C', 'B', 'A', 'B', 'C'],
    # 可以包含更多的协变量...
})

# 对分类变量进行独热编码
data = pd.get_dummies(data, columns=['癌症类型'], drop_first=True)

# 创建CoxPHFitter对象
cph = CoxPHFitter()

# 输入生存时间、事件发生和协变量数据
cph.fit(data, duration_col='生存时间', event_col='事件发生', formula='年龄 + 性别')

# 打印模型的参数
print(cph.summary)