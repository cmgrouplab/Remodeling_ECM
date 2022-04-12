import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation


parser = argparse.ArgumentParser()
parser.add_argument('--data-dir', type=str, default='.')
parser.add_argument('--save-dir', type=str, default='.')
parser.add_argument('--steps', type=int, default=10000,
                    help="number of movement steps")
parser.add_argument('--tick', type=int, default=10,
                    help="number of steps between each save file")
parser.add_argument('--radius', type=float, default=0.02,
                    help="radius of particles")
parser.add_argument('--length', type=float, default=1.0,
                    help="box length")
parser.add_argument('--num-grid', type=int, default=5,
                    help="num of grid along one side")
parser.add_argument('--interval', type=int, default=20,
                    help="time interval between animation frames in miliseconds")
parser.add_argument('--figsize', type=int, default=10,
                    help="figure size")
args = parser.parse_args()


def animate_data(patch_data, arrow_data, region, savedir):
    plt.rcParams['animation.html'] = 'html5'

    fig, ax = plt.subplots(figsize=(args.figsize, args.figsize))

    xmin, xmax, ymin, ymax = region
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect('equal')

    # def init():
    #     return []

    patches = []

    for position in patch_data[0]:
        patches.append(ax.add_patch(plt.Circle(
            (position[0], position[1]), args.radius, edgecolor='b', facecolor='y')))

    quivers = ax.quiver(FIBER_POSITIONS[:, 0], FIBER_POSITIONS[:, 1], arrow_data[0][:, 0], arrow_data[0][:, 1],
                        headlength=0, headwidth=1, pivot='mid', angles='xy', scale_units='xy', scale=1)

    def animate(i):
        for circle, position in zip(patches, patch_data[i]):
            circle.center = position[0], position[1]

        quivers.set_UVC(arrow_data[i][:, 0], arrow_data[i][:, 1])

        ax.set_title(f"Step={i*args.tick}")
        return patches, quivers

    anim = animation.FuncAnimation(fig, animate,
                                   frames=len(patch_data), interval=args.interval, blit=False)

    if not os.path.exists(savedir):
        os.mkdir(savedir)

    anim.save(os.path.join(savedir, 'Brownian.gif'),
              dpi=80, writer='imagemagick')


particle_files = [f"position{i}.csv" for i in range(
    0, args.steps+1, args.tick)]
particle_data = [np.loadtxt(os.path.join(args.data_dir, f), delimiter=',')
                 for f in particle_files]

GRID_SIZE = args.length / args.num_grid
FIBER_POSITIONS = np.array([[GRID_SIZE * (i+0.5), GRID_SIZE * (j+0.5)]
                            for i in range(args.num_grid) for j in range(args.num_grid)])

fiber_files = [f"fiberPosition{i}.csv" for i in range(
    0, args.steps + 1, args.tick)]
fiber_data = [np.loadtxt(os.path.join(args.data_dir, f), delimiter=',') * GRID_SIZE
              for f in fiber_files]

assert len(fiber_data[0]) == len(
    FIBER_POSITIONS), "Num of grids in data file not matching num_grid^2"

animate_data(particle_data, fiber_data, [
             0, args.length, 0, args.length], args.save_dir)
