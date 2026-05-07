import matplotlib.pyplot as plt
import numpy as np
import os

# =========================================================
# 📊 毕设第三章：PCG 求解器收敛历史分析图 (Semi-log Plot)
# =========================================================

# 智能寻路：自动在根目录和 build 目录中寻找数据文件
possible_paths = ["pcg_residual.txt", "build/pcg_residual.txt", "tests/pcg_residual.txt"]
file_path = None

for p in possible_paths:
    if os.path.exists(p):
        file_path = p
        break

if file_path is None:
    print("❌ 找不到 pcg_residual.txt 文件！")
    print("请确认你已经在 Docker 终端里执行了 ./tests/test_benchmark_iter 哇！")
    exit()

print(f"✅ 成功找到收敛数据文件: {file_path}")

# 读取两列数据：第一列是迭代次数，第二列是残差
try:
    data = np.loadtxt(file_path, skiprows=1)
    iters = data[:, 0]
    residuals = data[:, 1]
except Exception as e:
    print("数据读取失败，请检查 txt 文件格式！", e)
    exit()

# 开始绘制学术图表
fig, ax = plt.subplots(figsize=(8, 6))

# 核心：使用 semilogy (半对数Y轴) 绘制下降曲线
ax.semilogy(iters, residuals, color='#0033a0', linewidth=2.5, marker='o', markersize=3, label='PCG Relative Residual')

# 画一条红色的虚线，代表收敛容差 (1e-6)
tol = 1e-6
ax.axhline(y=tol, color='red', linestyle='--', linewidth=1.5, label=f'Convergence Tolerance ({tol})')

# 设置图表格式
ax.set_xlabel('Iteration Number', fontsize=12, fontweight='bold')
ax.set_ylabel('Relative Residual $||Ax-b|| / ||b||$ (Log Scale)', fontsize=12, fontweight='bold')
ax.set_title('Convergence History of Preconditioned Conjugate Gradient', fontsize=14, fontweight='bold')

# 开启精细网格辅助线
ax.grid(True, which="major", ls="-", alpha=0.5)
ax.grid(True, which="minor", ls=":", alpha=0.3)
ax.legend(fontsize=11)

plt.tight_layout()
plt.show()