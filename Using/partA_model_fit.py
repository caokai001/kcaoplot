#! bin/python
## 加载库
import pandas as pd
from fitter import Fitter
import matplotlib.pyplot as plt
from plotnine import *
from plotnine.data import *

## 读取文件，转成字典
data = {}
path=r"C:\Users\16926\Desktop\tf_all_background"
# with open(r"/public/home/kcao/Desktop/2018_NC_GS/model_select/PartA_background/tf_all_background") as f:
with open(path) as f:
    for line in f.readlines():
        line = line.strip().split()
        # print(line[0])
        # print(line[1:])
        data[str(line[0])] = line[1:]

## 字典转换成pandas
c = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in data.items()]))
c.head()


## fitter拟合,保存结果到test 字典对象中
model_fitting_dict = {}
for i in range(c.shape[1]):
    tf_filter_data = c.iloc[:, i].dropna().to_list()  # 去除NAN
    tf_name = c.columns[i]

    filter_data = list(map(eval,tf_filter_data))      # 字符串转换成数字
    f = Fitter(filter_data, distributions=['gamma', 'norm', 'uniform', "expon"])
    f.fit()
    # f.summary()
    model_fitting_dict[tf_name]=f._fitted_errors
    # dict 转成pd.dataFrame
model_fitting_pd = pd.DataFrame.from_dict(model_fitting_dict).T


## 画boxplot
plt.show(block=False)
# model_fitting_pd.boxplot()
# 修改数据绘图的格式
data_plot =pd.melt(model_fitting_pd,var_name="model",value_name="value")

base_plot= (ggplot(data_plot,aes('model','value'))+geom_boxplot(aes(fill = 'model'),notch = False)+coord_flip()+theme(legend_position="none"))
print(base_plot)
save_as_pdf_pages([base_plot + theme(figure_size=(8, 6))], filename ="model_select.pdf",path=r"C:\Users\16926\Desktop")



# test = {"A" : {'expon': 0.9344681445197056,
#  'gamma': 0.7468316227180921,
#  'norm': 0.7479262380924326,
#  'uniform': 0.8846573215044258},"B":{'expon': 0.9344681445197056,
#  'gamma': 0.7468316227180921,
#  'norm': 0.7479262380924326,
#  'uniform': 0.8846573215044258}}
