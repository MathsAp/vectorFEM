import sys
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

def read_vertices(filename):
    vertices = []
    with open(filename, 'r', encoding='utf-8') as f:
        for line in f:
            if line.strip() == '':
                continue
            x, y = map(float, line.split())
            vertices.append((x, y))
    return vertices

def read_triangles(filename):
    triangles = []
    with open(filename, 'r', encoding='utf-8') as f:
        for line in f:
            if line.strip() == '':
                continue
            i1, i2, i3 = map(int, line.split())
            triangles.append((i1, i2, i3))
    return triangles

def main():
    if len(sys.argv) != 3:
        print("❌ Неверное количество аргументов.")
        print("Использование: python draw_mesh.py vertices.txt triangles.txt")
        input("\nНажмите Enter, чтобы выйти...")
        sys.exit(1)

    vertices_file = sys.argv[1]
    triangles_file = sys.argv[2]

    vertices = read_vertices(vertices_file)
    triangles = read_triangles(triangles_file)

    fig, ax = plt.subplots()
    ax.set_aspect('equal')

    # Рисуем треугольники
    for tri in triangles:
        points = [vertices[i] for i in tri]
        polygon = Polygon(points, closed=True, edgecolor='black', facecolor='lightblue', alpha=0.7)
        ax.add_patch(polygon)

    # Отмечаем вершины и их номера
    for i, (x, y) in enumerate(vertices):
        ax.plot(x, y, 'ro', markersize=4)
        ax.text(x + 0.02, y + 0.02, str(i), fontsize=8, color='red')

    ax.autoscale_view()
    ax.grid(True)
    plt.title("Mesh")
    plt.show()

if __name__ == "__main__":
    main()
