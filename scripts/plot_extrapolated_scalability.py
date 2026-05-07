import matplotlib.pyplot as plt
import numpy as np

# 设置纯净学术字体
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif']
plt.rcParams['axes.unicode_minus'] = False    

# ==============================================================
# 👑 包含 O(N^2.15) 理论拟合外推的终极数据 👑
# ==============================================================
dof = [31800, 222480, 422400, 588120, 1121400]
time_pardiso = [0.612, 38.416, 140.658, 286.910, 1161.945]  # 112万使用 O(N^2.15) 严格拟合
time_cg      = [3.100, 60.837, 142.119, 158.994, 269.842]  
time_pcg     = [2.170, 86.623, 190.299, 297.385, 195.205]  

fig, ax = plt.subplots(figsize=(11, 7.5), dpi=150)

# 绘制曲线
ax.plot(dof, time_pardiso, marker='o', markersize=9, linewidth=2.5, color='#E15759', label='Intel MKL Pardiso (Theoretical Projection)')
ax.plot(dof, time_cg, marker='s', markersize=9, linewidth=2.5, color='#4E79A7', label='Standard CG (No Preconditioner)')
ax.plot(dof, time_pcg, marker='^', markersize=9, linewidth=2.5, linestyle='--', color='#76B7B2', label='Jacobi PCG (Diagonal Scaled)')

# 特殊标注 1.12M 处的 Pardiso (空心圆表示拟合，并保留 OOM 警告)
ax.plot(1121400, 1161.945, marker='o', markersize=14, markerfacecolor='white', markeredgecolor='#E15759', markeredgewidth=2)
ax.annotate('Theoretical: ~1162s\n(Actual: OOM Crash)', xy=(1121400, 1161.945), xytext=(800000, 1050),
            arrowprops=dict(facecolor='#E15759', shrink=0.05, width=1.5, headwidth=6),
            fontsize=12, color='#E15759', fontweight='bold', ha='center')

# 突出展示 PCG 在百万级的恐怖碾压优势
ax.annotate('PCG is ~6x Faster!\n(195.2s vs 1162s)', xy=(1121400, 195.205), xytext=(850000, 350),
            arrowprops=dict(facecolor='#76B7B2', shrink=0.05, width=1.5, headwidth=6),
            fontsize=12, color='#76B7B2', fontweight='bold', bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#76B7B2", alpha=0.9))

# 标注其他重要节点
ax.annotate('Golden Crossover\n(Direct vs CG Tie)', xy=(422400, 142.119), xytext=(150000, 300),
            arrowprops=dict(facecolor='#B8860B', shrink=0.05, width=1.5, headwidth=6),
            fontsize=11, color='#B8860B', fontweight='bold')

# 设置交叉区指示
ax.axvline(x=422400, color='#F2CF5B', linestyle='--', linewidth=2, alpha=0.6)

# 设置坐标轴与外观 (Y轴大幅拉高以容纳指数爆炸)
ax.set_xlabel('Total Degrees of Freedom (DOF)', fontsize=14, fontweight='bold')
ax.set_ylabel('Wall-clock Time (Seconds)', fontsize=14, fontweight='bold')
ax.set_title('Ultimate Solver Scalability: Exponential Explosion vs Linear Efficiency', fontsize=16, fontweight='bold', pad=20)

xticks = np.arange(0, 1300000, 200000)
ax.set_xticks(xticks)
ax.set_xticklabels([f"{int(x/1000)}k" for x in xticks], fontsize=12)
ax.set_ylim(0, 1250)
ax.tick_params(axis='both', labelsize=12)
ax.legend(fontsize=12, loc='upper left')
ax.grid(True, linestyle='--', alpha=0.5)

plt.tight_layout()
plt.savefig('theoretical_projection_scalability_curve.png', dpi=300)
print("✅ 带有理论外推的视觉震撼版图表已生成！请查看 'theoretical_projection_scalability_curve.png'！")
