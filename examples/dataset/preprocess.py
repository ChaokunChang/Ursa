import argparse


def parse_args():
    parser = argparse.ArgumentParser(description='Preprocess data..')
    parser.add_argument('--data-path')
    parser.add_argument('--save-path')
    parser.add_argument('--directed', action='store_true')

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_args()
    with open(args.data_path, 'r') as f:
        lines = f.readlines()

    with open(args.save_path, 'w') as f:
        graph = dict()
        for line in lines[1:]:
            line = line.strip().split()
            if line[0] == 'v':
                vid, vlabel = line[1:]
                graph[vid] = [vid, vlabel]
            else:
                vid_from, vid_to, elabel = line[1:]
                graph[vid_from] += [vid_to, elabel]
                if not args.directed:
                    graph[vid_to] += [vid_from, elabel]

        for k, v in graph.items():
            num_neighbors = len(v[2:]) // 2
            v = v[:2] + [str(num_neighbors)] + v[2:]
            f.write(' '.join(v) + '\n')
