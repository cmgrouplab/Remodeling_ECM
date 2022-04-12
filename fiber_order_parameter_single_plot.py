import os
import glob
import argparse
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg

parser = argparse.ArgumentParser()
parser.add_argument('--data-dir', type=str, default='.')
parser.add_argument('--save-dir', type=str, default='.')
parser.add_argument('--steps', type=int, default=10000,
                    help="number of movement steps")
parser.add_argument('--tick', type=int, default=100,
                    help="number of steps between each save file")
args = parser.parse_args()


def main(data_dir):
    files = [f"fiberPosition{i}.csv" for i in range(
        0, args.steps + 1, args.tick)]

    y = []
    directors = []
    for ff in files:
        file_label = ff[8:0-4]
        f = os.path.join(data_dir, ff)

        fiber = np.loadtxt(f, delimiter=',')[:, :2]
        # fiber (N, 2)
        # np.expand_dim(fiber, 2)  (N, 2, 1)
        # np.expand_dim(fiber, 1)  (N, 1, 2)
        A = np.mean(np.expand_dims(fiber, axis=2) *
                    np.expand_dims(fiber, axis=1), axis=0) * 1.5 - 0.5 * np.eye(2)
        ws, vs = scipy.linalg.eig(A)
        directors.append(vs[np.argmax(ws)])
        y.append(np.real(np.max(ws)))
    x = np.arange(0, args.steps + 1, args.tick)
    y = np.array(y)
    directors = np.array(directors)

    np.savetxt(os.path.join(data_dir, 'fiber_order_parameter.txt'),
               np.concatenate([x[:, np.newaxis], y[:, np.newaxis], directors], axis=1))
    fig, axes = plt.subplots(figsize=(10, 10))
    axes.plot(x, y)
    #axes.fill_between(x, y-y_err, y+y_err, alpha=0.4)
    # axes.set_xticks(bins)
    # axes.set_xticklabels(bins, rotation=45)
    axes.set_title(f'Fiber order parameter_{data_dir}')
    axes.set_xlabel('Step')
    axes.set_ylabel('order parameter')
    plt.ylim(0, 1)

    # axes.legend()
    plt.savefig(os.path.join(
        args.data_dir, f'fiber_order_parameter{"_".join(data_dir.split("/"))}.png'))


data_dirs = sorted(glob.glob(os.path.join(args.data_dir, 'ChangeCoeff*')))
print(data_dirs)


for data_dir in data_dirs:
    main(data_dir)
