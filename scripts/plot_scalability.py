import matplotlib.pyplot as plt
import numpy as np

# 设置纯净学术字体
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif']
plt.rcParams['axes.unicode_minus'] = False    

# ==============================================================
# 🟩 绝对真实的终端实测数据 (Ground Truth) 🟩
# ==============================================================
dof = [31800, 222480, 422400, 588120]
time_pardiso = [0.612, 38.416, 140.658, 286.910]   
time_cg      = [3.100, 60.837, 142.119, 158.994]  
time_pcg     = [2.170, 86.623, 190.299, 297.385]  

fig, ax = plt.subplots(figsize=(11, 7.5), dpi=150)

# 绘制曲线
ax.plot(dof, time_pardiso, marker='o', markersize=9, linewidth=2.5, color='#E15759', label='Intel MKL Pardiso (Direct)')
ax.plot(dof, time_cg, marker='s', markersize=9, linewidth=2.5, color='#4E79A7', label='Standard CG (No Preconditioner)')
ax.plot(dof, time_pcg, marker='^', markersize=9, linewidth=2.5, linestyle='--', color='#76B7B2', label='Jacobi PCG (Mathematical Convergence but High HW Overhead)')

# 标注 Pardiso 时间
for i, txt in enumerate(time_pardiso):
    offset_y = 15
    ax.annotate(f"{txt:.1f}s", (dof[i], time_pardiso[i]), textcoords="offset points", xytext=(-10, offset_y), ha='right', fontsize=11, color='#E15759', fontweight='bold')

# 标注 CG 时间 (高亮与Pardiso的死斗)
for i, txt in enumerate(time_cg):
    offset_y = -18
    ax.annotate(f"{txt:.1f}s", (dof[i], time_cg[i]), textcoords="offset points", xytext=(0, offset_y), ha='center', fontsize=11, color='#4E79A7', fontweight='bold')

# 标注 PCG 时间
for i, txt in enumerate(time_pcg):
    offset_y = 15
    ax.annotate(f"{txt:.1f}s", (dof[i], time_pcg[i]), textcoords="offset points", xytext=(12, offset_y), ha='left', fontsize=11, color='#76B7B2', fontweight='bold')

# 标注核心区域：性能交叉点！
ax.axvspan(400000, 440000, color='#F2CF5B', alpha=0.2, label='Scalability Crossover Zone')
ax.annotate('Golden Crossover\n(Direct vs CG Tie)', xy=(422400, 141), xytext=(280000, 220),
            arrowprops=dict(facecolor='#B8860B', shrink=0.05, width=1.5, headwidth=6),
            fontsize=12, color='#B8860B', fontweight='bold')

# 设置坐标轴与外观
ax.set_xlabel('Total Degrees of Freedom (DOF)', fontsize=14, fontweight='bold')
ax.set_ylabel('Wall-clock Time (Seconds)', fontsize=14, fontweight='bold')
ax.set_title('Solver Scalability: Direct vs. Iterative Benchmark (Ground Truth Data)', fontsize=16, fontweight='bold', pad=20)

xticks = np.arange(0, 700000, 100000)
ax.set_xticks(xticks)
ax.set_xticklabels([f"{int(x/1000)}k" for x in xticks], fontsize=12)
ax.tick_params(axis='both', labelsize=12)
ax.legend(fontsize=11, loc='upper left')
ax.grid(True, linestyle='--', alpha=0.5)

plt.tight_layout()
plt.savefig('real_data_scalability_curve_ground_truth.png', dpi=300)
print("✅ 终极严谨的真实数据点折线图已生成！请在 Windows 文件夹中查看！")
