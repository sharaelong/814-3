import matplotlib.pyplot as plt

def visualize_clusters(filename):
    with open(filename, 'r') as file:
        n, k = map(int, file.readline().split())
        points = [list(map(int, line.split())) for line in file.readlines()[:n]]
        clusters = [list(map(int, line.split()[1:])) for line in file.readlines()]

    # Extract the coordinates of the points
    x = [point[0] for point in points]
    y = [point[1] for point in points]

    # Plot the points
    plt.scatter(x, y, color='black')

    # Assign a different color to each cluster
    colors = ['red', 'green', 'blue', 'orange', 'purple', 'yellow']
    for i in range(k):
        cluster_points = [points[idx - 1] for idx in clusters[i]]
        cluster_x = [point[0] for point in cluster_points]
        cluster_y = [point[1] for point in cluster_points]
        plt.scatter(cluster_x, cluster_y, color=colors[i % len(colors)])

    # Set plot title and labels
    plt.title('Point Clusters Visualization')
    plt.xlabel('X')
    plt.ylabel('Y')

    # Show the plot
    plt.show()

# Usage example
filename = 'data.in'
visualize_clusters(filename)
