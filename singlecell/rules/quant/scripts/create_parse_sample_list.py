import yaml
import sys

if __name__=="__main__":
    conf_fn = sys.argv[1]
    out = sys.argv[2]

    with open(conf_fn) as fh:
        c = yaml.safe_load(fh)
    with open(out, 'w+') as fh:
        for k, v in c['wells'].items():
            fh.write("{} {}\n".format(k, v['Wells']))
