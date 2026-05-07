import matplotlib.pyplot as plt
import numpy as np

# 设置纯净学术字体
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif']
plt.rcParams['axes.unicode_minus'] = False    

# ==============================================================
# 🚨 100% 纯正带孔板真实数据 (已彻底纠正图1的混淆)
# ==============================================================
# X轴：系统总自由度 (DOF)
dof = [31800, 222480, 588120]

# Y轴：求解耗时 (Seconds)
time_pardiso = [0.612, 38.416, 286.910]   
time_cg = [3.100, 60.837, 158.994]        
time_pcg = [2.170, 86.623, 297.385]       

# ==============================================================

fig, ax = plt.subplots(figsize=(10.5, 7), dpi=150)

# 绘制折线
ax.plot(dof, time_pardiso, marker='o', markersize=9, linewidth=2.5, color='#E15759', label='Intel MKL Pardiso (Direct)')
ax.plot(dof, time_cg, marker='s', markersize=9, linewidth=2.5, color='#4E79A7', label='Standard CG (Iterative)')
ax.plot(dof, time_pcg, marker='^', markersize=9, linewidth=2.5, linestyle='--', color='#76B7B2', label='Jacobi PCG (Iterative)')

# 标注真实时间
for i, txt in enumerate(time_pardiso):
    offset_y = 12 if i == 2 else 15
    ax.annotate(f"{txt:.1f}s", (dof[i], time_pardiso[i]), textcoords="offset points", xytext=(-10, offset_y), ha='right', fontsize=11, color='#E15759', fontweight='bold')

for i, txt in enumerate(time_cg):
    ax.annotate(f"{txt:.1f}s", (dof[i], time_cg[i]), textcoords="offset points", xytext=(12, -15), ha='left', fontsize=11, color='#4E79A7', fontweight='bold')

for i, txt in enumerate(time_pcg):
    offset_y = -20 if i == 2 else 12
    ax.annotate(f"{txt:.1f}s", (dof[i], time_pcg[i]), textcoords="offset points", xytext=(12, offset_y), ha='left', fontsize=11, color='#76B7B2', fontweight='bold')

# 绘制交叉点提示区 (Crossover Zone)
# Pardiso 和 CG 的真实交叉点发生在 22万 到 58万 自由度之间
ax.axvspan(300000, 480000, color='#F2CF5B', alpha=0.15, label='Performance Crossover Zone')
ax.text(390000, 100, "Crossover\nPoint\n(O(N) beats O(N²))", ha='center', va='center', fontsize=12, fontweight='bold', color='#B8860B')

# 设置坐标轴
ax.set_xlabel('Total Degrees of Freedom (DOF)', fontsize=14, fontweight='bold')
ax.set_ylabel('Wall-clock Time (Seconds)', fontsize=14, fontweight='bold')
ax.set_title('Solver Scalability: Direct vs. Iterative Methods (Strict Control: Plate with Hole)', fontsize=16, fontweight='bold', pad=20)

# 优化X轴刻度显示格式
xticks = np.arange(0, 700000, 100000)
ax.set_xticks(xticks)
ax.set_xticklabels([f"{int(x/1000)}k" for x in xticks], fontsize=12)
ax.tick_params(axis='both', labelsize=12)

ax.legend(fontsize=12, loc='upper left')
ax.grid(True, linestyle='--', alpha=0.5)

plt.tight_layout()
plt.savefig('real_data_scalability_curve_corrected.png', dpi=300)
print("✅ 纯正带孔板折线图已完美纠正并生成！请查看 real_data_scalability_curve_corrected.png")
