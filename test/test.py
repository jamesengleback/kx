import sys
import kx

def main(args):
    targets = args 
    if targets == []:
        targets = ['1.1.1.1']
    for i in targets:
        print(kx.kx(i, 'data'))

if __name__ == '__main__':
    main(sys.argv[1:])
