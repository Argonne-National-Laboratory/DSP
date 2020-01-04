from coinor.grumpy import BBTree

bt = BBTree()
file = open('vbc.out', 'r')
for line in file:
    bt.ProcessLine(line)

bt.set_display_mode('file')
bt.display()
