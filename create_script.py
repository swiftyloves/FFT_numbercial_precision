

inp = input('Enter the highest power of x: ')
f = open('coefficient_polynomial' +inp+'.txt','w')

power = int(input('Enter the power of 10 of coefficient: '));



for i in range(int(inp)):
    f.write(inp + "\n")
    for j in range(int(inp)):
        if i == j:
          f.write(str(pow(10, power)) + " ")
        else:
          f.write("0 ")
    f.write("\n")
