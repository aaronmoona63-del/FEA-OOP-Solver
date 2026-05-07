import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# 设置全局学术风格
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 12

# 增加画布比例，防止边缘裁剪
fig = plt.figure(figsize=(12, 9))
ax = fig.add_subplot(111, projection='3d')

# 几何尺寸 L=50, W=10, H=10
L, W, H = 50, 10, 10

# 定义顶点坐标
vertices = np.array([
    [0, 0, 0], [L, 0, 0], [L, W, 0], [0, W, 0],
    [0, 0, H], [L, 0, H], [L, W, H], [0, W, H]
])

# 定义六个面
faces = [
    [vertices[0], vertices[1], vertices[5], vertices[4]], # Front
    [vertices[1], vertices[2], vertices[6], vertices[5]], # Right (Load)
    [vertices[2], vertices[3], vertices[7], vertices[6]], # Back
    [vertices[3], vertices[0], vertices[4], vertices[7]], # Left (Fixed)
    [vertices[0], vertices[1], vertices[2], vertices[3]], # Bottom
    [vertices[4], vertices[5], vertices[6], vertices[7]]  # Top
]

# 1. 绘制主体长方体
poly = Poly3DCollection(faces, facecolors='#adcbe3', linewidths=1.2, edgecolors='black', alpha=0.5)
ax.add_collection3d(poly)

# 2. 绘制红色节点
ax.scatter(vertices[:, 0], vertices[:, 1], vertices[:, 2], color='red', s=40, edgecolors='black', zorder=10)

# 3. 绘制绿色固定面 (X=0 平面)
fixed_face_verts = [[vertices[3], vertices[0], vertices[4], vertices[7]]]
fixed_poly = Poly3DCollection(fixed_face_verts, facecolors='green', alpha=0.3, hatch='////', edgecolor='green')
ax.add_collection3d(fixed_poly)

# 在固定面上增加三角形约束标记 (3x3分布)
fy, fz = np.meshgrid(np.linspace(0, W, 3), np.linspace(0, H, 3))
fx = np.zeros_like(fy)
ax.scatter(fx, fy, fz, color='green', marker='>', s=100, alpha=0.9, zorder=5)

# 4. 绘制橙色拉伸载荷箭头 (X=L 平面)
ly, lz = np.meshgrid(np.linspace(0, W, 3), np.linspace(0, H, 3))
lx = np.full_like(ly, L)
ax.quiver(lx, ly, lz, 15, 0, 0, color='#ff7f0e', linewidth=2.5, arrow_length_ratio=0.3)

# 5. 尺寸标注 (优化引出线距离以防裁剪)
# 标注 L
ax.plot([0, L], [-8, -8], [0, 0], color='black', lw=1.2, marker='|')
ax.text(L/2, -13, 0, '$L = 50$ mm', ha='center', style='italic')

# 标注 W (将其放置在远离中心的位置)
ax.plot([L+2, L+2], [0, W], [H+8, H+8], color='black', lw=1.2, marker='|')
ax.text(L+4, W/2, H+10, '$W = 10$ mm', ha='left', va='center', style='italic')

# 标注 H
ax.plot([-8, -8], [W, W], [0, H], color='black', lw=1.2, marker='_')
ax.text(-14, W, H/2, '$H = 10$ mm', ha='right', va='center', rotation=90, style='italic')

# 6. 自定义坐标轴 (移至左侧空白区)
off = -20
ax.quiver(off, off, 0, 15, 0, 0, color='black', lw=1.5, arrow_length_ratio=0.2)
ax.text(off+18, off, 0, 'X', weight='bold')
ax.quiver(off, off, 0, 0, 15, 0, color='black', lw=1.5, arrow_length_ratio=0.2)
ax.text(off, off+18, 0, 'Y', weight='bold')
ax.quiver(off, off, 0, 0, 0, 15, color='black', lw=1.5, arrow_length_ratio=0.2)
ax.text(off, off, 18, 'Z', weight='bold')

# 7. 文字标签 (增加偏移量，防止被切)
ax.text(L+20, W/2, -5, 'Traction: 100 MPa', color='#ff7f0e', fontsize=11, weight='bold', ha='center')
ax.text(0, W/2, H+18, 'Fixed Surface ($Z=0$): $U_x=U_y=U_z=0$', color='green', fontsize=11, weight='bold', ha='center')

# 8. 视角与显示范围控制 (关键：扩大边界防止裁剪)
ax.set_box_aspect((50, 20, 20))
ax.set_axis_off()

# 扩大坐标轴极限，为标注留出空间
ax.set_xlim([-25, 75])
ax.set_ylim([-20, 30])
ax.set_zlim([-10, 35])

# 视角微调
ax.view_init(elev=25, azim=-65)

# 强制调整子图边距
plt.subplots_adjust(left=0, right=1, bottom=0, top=1)

# 导出图片
output_path = 'FEA_Model_Complete.png'
plt.savefig(output_path, dpi=300, bbox_inches='tight', pad_inches=0.5)
plt.show()

print(f"✅ 图片已成功保存至: {output_path}，所有标注均已完整保留！")