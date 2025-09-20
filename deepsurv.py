import deepsurv

import sys
sys.path.append('E:/article4.5/smrsp2/DeepSurv/deepsurv')
import deepsurv
from deepsurv import deep_surv
from deepsurv_logger import DeepSurvLogger, TensorboardLogger
import utils
import viz
import numpy as np
import pandas as pd
import lasagne
import matplotlib
import matplotlib.pyplot as plt
import theano
# 读入example_data.csv
train_dataset_fp = 'E:/article4.5/smrsp2/train_data.csv'
train_df = pd.read_csv(train_dataset_fp)
train_df.head()
test_dataset_fp = 'E:/article4.5/smrsp2/vali_data.csv'
test_df = pd.read_csv(test_dataset_fp)
test_df.head()


def dataframe_to_deepsurv_ds(df, event_col = 'Event', time_col = 'Time'):
    e = df[event_col].values.astype(np.int32)
    t = df[time_col].values.astype(np.float32)
    x_df = df.drop([event_col, time_col], axis = 1)
    x = x_df.values.astype(np.float32)
    return {
        'x' : x,
        'e' : e,
        't' : t
    }

# event_col、time_col分别定义为自己数据的表头（在这里分别是Event和Time）
# 使用该函数对训练数据集进行转换，同样可以用于测试数据集、验证数据集
train_data = dataframe_to_deepsurv_ds(train_df, event_col = 'Event', time_col= 'Time')
test_data = dataframe_to_deepsurv_ds(test_df, event_col = 'Event', time_col= 'Time')


hyperparams = {
    'L2_reg': 10.0,
    'batch_norm': True,
    'dropout': 0.4,
    'hidden_layers_sizes': [25, 25],
    'learning_rate': 1e-05,
    'lr_decay': 0.001,
    'momentum': 0.9,
    'n_in': train_data['x'].shape[1],
    'standardize': True
}



#创建一个实例
model = deep_surv.DeepSurv(**hyperparams)

#（可选） 利用 TensorBoard 来监控训练和验证
# 这个步骤是可选
# 如果你不想使用tensorboard logger，请取消下面一行的注释"#"，并注释掉其他三行
# logger = None
experiment_name = 'test_experiment_sebastian'
logdir = './logs/tensorboard/'
logger = TensorboardLogger(experiment_name, logdir=logdir)

# 开始训练我们的模型
update_fn=lasagne.updates.nesterov_momentum # 你目前使用的optimizer，可以检查以下网址获得更多optimizer：http://lasagne.readthedocs.io/en/latest/modules/updates.html
# 定义epoch迭代数量，适量即可
n_epochs =1000

# 如果有验证数据valid_data，可以将其作为第二个参数添加到下方函数中：
# model.train(train_data, valid_data, n_epochs=500)
metrics = model.train(train_data, n_epochs=n_epochs, logger=logger, update_fn=update_fn)

import torch

torch.save(model, 'E:/article4.5/smrsp2/model')
model2=torch.load('E:/article4.5/smrsp2/model')


deepsurv.plot_log(metrics)

risk_train=model.predict_risk(train_data["x"])
risk_train=pd.DataFrame(risk_train)

risk_train.to_csv('E:/article4.5/smrsp2/risk_train.csv')


risk_test=model.predict_risk(test_data["x"])
risk_test=pd.DataFrame(risk_test)

risk_test.to_csv('E:/article4.5/smrsp2/risk_vali.csv')

train_df = pd.read_csv(train_dataset_fp)


