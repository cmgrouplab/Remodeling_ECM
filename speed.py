
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

parser = argparse.ArgumentParser()
parser.add_argument('--data-dir', type=str, default='.')
parser.add_argument('--steps', type=int, default=30000,
                    help="number of movement steps")
parser.add_argument('--tick', type=int, default=1000,
                    help="number of steps between each save file")
parser.add_argument('--delta_r', type=float, default=0.002,
                    help="number of steps between each save file")
args = parser.parse_args()


def speed(data, delta_r):
    positions = data[:, :2]
    velocities = data[:, 2:]
    aver_speed = []
    center = np.array([0.5, 0.5])
    for r in np.arange(0.05, 0.51, 2 * delta_r):
        total_speed = []
        for i in range(len(positions)):
            distance = np.linalg.norm(positions[i]-center)
            if distance < r + delta_r and distance >= r - delta_r:
                s = np.linalg.norm(velocities[i])
                total_speed.append(s)

        aver_speed.append([r, np.mean(total_speed)])
    return aver_speed


def animate_data(data, savedir):
    plt.rcParams['animation.html'] = 'html5'

    fig, ax = plt.subplots(figsize=(10, 10))

    ax.set_xlim(0, 0.51)
    ax.set_ylim(0, 5e-6)
    # ax.set_aspect('equal')
    line, = ax.plot([], [])

    def init():
        return line,

    def animate(i):
        x = data[i][:, 0]
        y = data[i][:, 1]
        line.set_data(x, y)

        ax.set_title(f"Step={i*args.tick}")
        return line,

    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=len(data), interval=args.interval, blit=True)

    if not os.path.exists(savedir):
        os.mkdir(savedir)

    anim.save(os.path.join(savedir, 'speed_plot.gif'),
              dpi=80, writer='imagemagick')


files = [f"position{i}.csv" for i in range(0, args.steps+1, args.tick)]
all_data = [np.loadtxt(os.path.join(args.data_dir, f),
                       delimiter=',') for f in files]
all_speeds = [speed(data, args.delta_r) for data in all_data]

animate_data(all_speeds, args.save_dir)
