import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif']
plt.rcParams['axes.unicode_minus'] = False

# 真实数据
scales = ['Small (31k)', 'Medium (222k)', 'Large (588k)']
pardiso = [0.612, 38.416, 286.910]
cg = [3.100, 60.837, 158.994]
pcg = [2.170, 86.623, 297.385]

x = np.arange(len(scales))
width = 0.25

fig, ax = plt.subplots(figsize=(11, 7), dpi=150)

rects1 = ax.bar(x - width, pardiso, width, label='Intel MKL Pardiso', color='#E15759', edgecolor='black', alpha=0.9)
rects2 = ax.bar(x, cg, width, label='Standard CG (Proposed)', color='#4E79A7', edgecolor='black', alpha=0.9)
rects3 = ax.bar(x + width, pcg, width, label='Jacobi PCG', color='#76B7B2', edgecolor='black', alpha=0.9)

ax.set_ylabel('Execution Time (Seconds)', fontsize=12, fontweight='bold')
ax.set_title('Performance Comparison Across Mesh Scales', fontsize=14, fontweight='bold', pad=15)
ax.set_xticks(x)
ax.set_xticklabels(scales, fontsize=11)
ax.legend()

# 标注数据
def autolabel(rects):
    for rect in rects:
        height = rect.get_height()
        ax.annotate(f'{height:.1f}s', xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=9)

autolabel(rects1); autolabel(rects2); autolabel(rects3)

# 增加对数坐标缩略图 (Inset) 处理小规模数据看不清的问题
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
ax_ins = inset_axes(ax, width="30%", height="30%", loc='upper center', borderpad=3)
ax_ins.bar(x - width, pardiso, width, color='#E15759', edgecolor='black')
ax_ins.bar(x, cg, width, color='#4E79A7', edgecolor='black')
ax_ins.bar(x + width, pcg, width, color='#76B7B2', edgecolor='black')
ax_ins.set_yscale('log')
ax_ins.set_title('Log Scale View', fontsize=9)
ax_ins.set_xticks(x)
ax_ins.set_xticklabels(['31k', '222k', '588k'], fontsize=8)

plt.tight_layout()
plt.savefig('final_performance_bar_chart.png')
print("✅ 终极真实数据柱状图已生成：final_performance_bar_chart.png")
