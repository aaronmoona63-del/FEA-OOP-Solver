import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib

# 设置学术字体风格 (防止中文乱码)
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.family'] = 'sans-serif'

# 1. 节点坐标 (基于 MATLAB 算例，乘以 1000 转换为 mm)
nodes = np.array([
    [1, 0, 0], [1, 0, 1], [0, 0, 1], [0, 0, 0],
    [1, 1, 1], [0, 2, 1], [0, 2, 0], [1, 1, 0],
    [2, 1, 1], [2, 1, 0], [2, 0, 1], [2, 0, 0]
]) * 1000.0

# 2. 单元拓扑连接 (将 MATLAB 的 1-based 转换为 Python 的 0-based)
# 注意：严格按照标准的 Hex8 逆时针节点顺序重排，以保证面片渲染正确
el1 = [7, 4, 5, 6, 0, 1, 2, 3] 
el2 = [1, 10, 11, 0, 4, 8, 9, 7] 

fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')

# 绘制带有半透明面片的 Hex8 单元
def draw_hex(ax, el, nodes, color='#3498db', edge_color='#2c3e50', alpha=0.15):
    # Hex8 的 6 个面
    faces = [
        [el[0], el[1], el[2], el[3]], # 底面
        [el[4], el[5], el[6], el[7]], # 顶面
        [el[0], el[1], el[5], el[4]], # 前面
        [el[2], el[3], el[7], el[6]], # 后面
        [el[1], el[2], el[6], el[5]], # 右面
        [el[0], el[3], el[7], el[4]]  # 左面
    ]
    # 绘制实体面片和加粗边缘
    poly3d = [[nodes[vert] for vert in face] for face in faces]
    ax.add_collection3d(Poly3DCollection(poly3d, facecolors=color, 
                                         linewidths=1.5, edgecolors=edge_color, alpha=alpha))

# 渲染两个单元 (分别用两种高级学术蓝/绿配色)
draw_hex(ax, el1, nodes, color='#2980b9')
draw_hex(ax, el2, nodes, color='#1abc9c')

# 3. 绘制红色的物理节点并标注节点号
ax.scatter(nodes[:,0], nodes[:,1], nodes[:,2], color='#c0392b', s=80, zorder=5, edgecolors='white', linewidth=1.5)
for i, p in enumerate(nodes):
    # 节点编号稍微向上偏移，防止重叠
    ax.text(p[0]+40, p[1]+40, p[2]+60, f'{i+1}', color='black', fontsize=13, fontweight='bold', zorder=6)

# 4. 绘制左侧固定边界条件 (X=0 平面上的节点: 3, 4, 6, 7 -> 对应 Python 索引 2, 3, 5, 6)
fixed_nodes = [2, 3, 5, 6]
for fn in fixed_nodes:
    p = nodes[fn]
    # 画绿色向左的三角形代表固支边界
    ax.scatter(p[0], p[1], p[2], color='#27ae60', marker='<', s=250, zorder=4)

ax.text(-900, 1000, 500, 'Fixed BC\n($U_x=U_y=U_z=0$)', color='#27ae60', 
        fontsize=15, fontweight='bold', ha='center')

# 5. 绘制右侧拉伸载荷 (X=2000 平面的受力节点: 11, 9, 10, 12 -> 对应 Python 索引 10, 8, 9, 11)
traction_nodes = [10, 8, 9, 11]
for tn in traction_nodes:
    p = nodes[tn]
    # 画橙色箭头表示拉伸载荷，方向严格指向 +X (u=800, v=0, w=0)
    ax.quiver(p[0], p[1], p[2], 800, 0, 0, color='#e67e22', 
              arrow_length_ratio=0.25, linewidth=3.5, zorder=5)

# ✨标明载荷大小
ax.text(2800, 500, 500, 'Traction Load\n$t_x = 100$ MPa', color='#e67e22', 
        fontsize=16, fontweight='bold')

# 6. 设置视角与坐标轴样式
# 🌟核心修复：强制锁定物理真实比例 (X:2000, Y:2000, Z:1000) 比例为 2:2:1
ax.set_box_aspect([2, 2, 1])

ax.set_xlim([-500, 3000])
ax.set_ylim([-500, 2500])
ax.set_zlim([-500, 1500])

ax.set_xlabel('X 坐标 (mm)', fontsize=13, labelpad=10)
ax.set_ylabel('Y 坐标 (mm)', fontsize=13, labelpad=10)
ax.set_zlabel('Z 坐标 (mm)', fontsize=13, labelpad=10)

# 去除灰色背景墙，纯白底色 SCI 顶级期刊质感
ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False
ax.xaxis.pane.set_edgecolor('white')
ax.yaxis.pane.set_edgecolor('white')
ax.zaxis.pane.set_edgecolor('white')
ax.grid(True, linestyle='-.', alpha=0.4, color='#7f8c8d')

# 绝佳观察视角
ax.view_init(elev=22, azim=-52)
plt.tight_layout()

# 导出高清图片
output_filename = 'Figure_3.8_Hex8_Benchmark.png'
plt.savefig(output_filename, dpi=600, bbox_inches='tight', transparent=True)
print(f"✅ 图片生成成功！已保存为: {output_filename}")