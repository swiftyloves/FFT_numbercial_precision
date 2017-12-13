

inp = input('Enter the highest power of x: ')
amount = input('Enter the amount of x: ')
f = open('coefficient_polynomial_all1_' +inp+'_'+amount+'.txt','w')

# power = int(input('Enter the power of 10 of coefficient: '));


f.write(inp + "\n")
for i in range(int(amount)):
    if i < int(amount):
        f.write("1 ")
    else:
        f.write("0 ")
        # if i == j:
        #   f.write(str(pow(10, power)) + " ")
        # else:
        #   f.write("0 ")
    # f.write("\n")
