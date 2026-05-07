import matplotlib.pyplot as plt
import numpy as np

# 纯净学术字体设置 (防报错版)
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif']
plt.rcParams['axes.unicode_minus'] = False    

# ==============================================================
# 🚨 你的所有真实心血数据！
# ==============================================================
# 第一组：小规模带孔板 (~31k DOF)
time_31k = [0.600, 3.100, 2.170]       # 耗时 [Pardiso, CG, PCG]
iter_31k = [0, 1646, 843]              # 迭代次数 [0, CG, PCG]

# 第二组：超大规模带孔板 (~588k DOF) 
time_588k = [286.910, 158.994, 297.385] # 刚刚跑出的耗时 [Pardiso, CG, PCG]
iter_588k = [0, 5193, 5527]             # 刚刚跑出的迭代 [0, CG, PCG]
# ==============================================================

labels = ['Intel MKL Pardiso\n(Direct Solver)', 'Standard CG\n(No Preconditioner)', 'Jacobi PCG\n(Proposed Method)']
x = np.arange(len(labels))
width = 0.35  

# 创建画布
fig, ax = plt.subplots(figsize=(11, 7), dpi=150)

# 绘制双柱状图
rects1 = ax.bar(x - width/2, time_31k, width, label='Small Scale (~31k DOF)', color='#98ABC5', edgecolor='black', linewidth=1.2)
rects2 = ax.bar(x + width/2, time_588k, width, label='Large Scale (~588k DOF)', color='#2B5B84', edgecolor='black', linewidth=1.2)

# 设置对数坐标系 (Log Scale)，这样 0.6秒 和 300秒 才能和谐地出现在同一张图里
ax.set_yscale('log')
ax.set_ylabel('Wall-clock Time (Seconds) [Log Scale]', fontsize=14)
ax.set_title('Solver Scalability Benchmark: Small vs. Large Scale Mesh', fontsize=16, pad=20)
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=12)
ax.legend(fontsize=12, loc='upper left')
ax.grid(axis='y', linestyle='--', alpha=0.5, which='both')
ax.tick_params(axis='y', labelsize=12)
ax.tick_params(axis='x', labelsize=12)

# 自动标注时间和迭代次数
def autolabel(rects, times, iters):
    for rect, t, it in zip(rects, times, iters):
        height = rect.get_height()
        # 标注耗时
        ax.text(rect.get_x() + rect.get_width()/2., height * 1.15,
                f'{t:.1f} s', ha='center', va='bottom', fontsize=11)
        
        # 标注迭代次数
        text = "Direct\n(No Iter)" if it == 0 else f"{it} Iters"
        ax.text(rect.get_x() + rect.get_width()/2., height * 0.5, text,
                ha='center', va='center', fontsize=10, color='white',
                bbox=dict(facecolor='black', alpha=0.5, edgecolor='none', pad=2))

autolabel(rects1, time_31k, iter_31k)
autolabel(rects2, time_588k, iter_588k)

# 留出顶部空间
ax.set_ylim(0.1, max(time_588k) * 5)

plt.tight_layout()
plt.savefig('ultimate_scalability_benchmark.png', dpi=300)
print("✅ 终极标度律双柱对比图已成功生成！请下载查看 ultimate_scalability_benchmark.png")
