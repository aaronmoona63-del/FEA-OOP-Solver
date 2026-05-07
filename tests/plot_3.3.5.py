import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib

# 设置学术字体风格
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman']

# 1. 节点坐标 (基于你的 MATLAB 算例，乘以 1000 转换为 mm)
nodes = np.array([
    [1, 0, 0], [1, 0, 1], [0, 0, 1], [0, 0, 0],
    [1, 1, 1], [0, 2, 1], [0, 2, 0], [1, 1, 0],
    [2, 1, 1], [2, 1, 0], [2, 0, 1], [2, 0, 0]
]) * 1000.0

# 2. 单元连接 (1-based 转换为 Python 的 0-based)
el1 = [8, 5, 6, 7, 1, 2, 3, 4]
el2 = [2, 11, 12, 1, 5, 9, 10, 8]
el1 = [x-1 for x in el1]
el2 = [x-1 for x in el2]

# 创建 3D 图形
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# 绘制 Hex8 单元网格线
def draw_hex(ax, el, nodes, color='#2c3e50', alpha=0.8):
    # Hex8 的 12 条拓扑边
    hex_edges = [
        (0,1), (1,2), (2,3), (3,0), # 底面
        (4,5), (5,6), (6,7), (7,4), # 顶面
        (0,4), (1,5), (2,6), (3,7)  # 柱面边
    ]
    for e in hex_edges:
        p1 = nodes[el[e[0]]]
        p2 = nodes[el[e[1]]]
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], 
                color=color, lw=2.0, alpha=alpha, zorder=3)

draw_hex(ax, el1, nodes)
draw_hex(ax, el2, nodes)

# 3. 绘制节点并标注节点号
ax.scatter(nodes[:,0], nodes[:,1], nodes[:,2], color='#e74c3c', s=50, zorder=5)
for i, p in enumerate(nodes):
    # 节点号稍微偏移一点，避免被点挡住
    ax.text(p[0]+50, p[1]+50, p[2]+50, f'{i+1}', 
            color='black', fontsize=12, fontweight='bold', zorder=6)

# 4. 绘制左侧固定边界条件 (X=0 平面的节点: 3, 4, 6, 7 -> 对应索引 2, 3, 5, 6)
fixed_nodes = [2, 3, 5, 6]
for fn in fixed_nodes:
    p = nodes[fn]
    # 使用向左的绿色三角形表示固支 (Fixed BC)
    ax.scatter(p[0], p[1], p[2], color='#27ae60', marker='<', s=180, zorder=4)

ax.text(-800, 1000, 500, 'Fixed BC\n(U = 0)', color='#27ae60', 
        fontsize=14, fontweight='bold', ha='center')

# 5. 绘制右侧拉伸载荷 (X=2000 平面的节点: 11, 9, 10, 12 -> 对应索引 10, 8, 9, 11)
traction_nodes = [10, 8, 9, 11]
for tn in traction_nodes:
    p = nodes[tn]
    # 画橙色箭头表示拉伸载荷
    ax.quiver(p[0], p[1], p[2], 600, 0, 0, color='#e67e22', 
              arrow_length_ratio=0.3, linewidth=2.5, zorder=4)

ax.text(2600, 500, 500, 'Traction', color='#e67e22', 
        fontsize=14, fontweight='bold')

# 6. 设置视角、坐标轴与样式
ax.set_xlabel('X (mm)', fontsize=12, labelpad=10)
ax.set_ylabel('Y (mm)', fontsize=12, labelpad=10)
ax.set_zlabel('Z (mm)', fontsize=12, labelpad=10)

ax.set_xlim([-400, 2800])
ax.set_ylim([-200, 2200])
ax.set_zlim([-200, 1200])

# 去除自带的灰色背景面板，使其符合SCI论文风格
ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False
ax.xaxis.pane.set_edgecolor('white')
ax.yaxis.pane.set_edgecolor('white')
ax.zaxis.pane.set_edgecolor('white')
ax.grid(True, linestyle='--', alpha=0.5)

# 调整一个能够完美展现两单元拓扑的视角
ax.view_init(elev=20, azim=-55)
plt.tight_layout()

# 保存高清插图
output_name = 'Figure_3.8_Hex8_Model.png'
plt.savefig(output_name, dpi=600, bbox_inches='tight', transparent=True)
print(f"✅ 图片已成功保存为: {output_name}")
# plt.show() # 取消注释可在窗口直接预览