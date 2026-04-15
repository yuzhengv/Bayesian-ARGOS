# %%
import numpy as np
import pysindy as ps
from scipy.integrate import solve_ivp

# %%
# 定义总人口数和其他参数
N = 10000
beta = 0.3  # 初始假设值
sigma = 1 / 5.2  # 潜伏率
gamma = 1 / 18  # 康复率

# 初始条件
I0 = 1
E0 = 10
R0 = 0
S0 = N - I0 - E0 - R0

# 时间范围
t = np.linspace(0, 100, 10000)


# 定义 SEIR 模型的微分方程
def seir_model(t, y):
    S, E, I, R = y
    dS_dt = -beta * S * I / N
    dE_dt = beta * S * I / N - sigma * E
    dI_dt = sigma * E - gamma * I
    dR_dt = gamma * I
    return [dS_dt, dE_dt, dI_dt, dR_dt]


# 求解微分方程
solution = solve_ivp(seir_model, [t[0], t[-1]], [S0, E0, I0, R0], t_eval=t)
S, E, I, R = solution.y

# %%
# 构建 PySINDy 模型
library = ps.PolynomialLibrary(degree=2)  # 使用多项式库
optimizer = ps.STLSQ(threshold=0.005)  # 使用稀疏回归优化器
model = ps.SINDy(feature_library=library, optimizer=optimizer)

# 准备数据
data = np.vstack((S, E, I, R)).T  # 将数据组合成一个矩阵
model.fit(data, t=t, multiple_trajectories=False)

# 打印发现的模型
model.print()

# %%
