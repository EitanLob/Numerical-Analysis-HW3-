import sympy

def bitsection(polynom, small_x, large_x, eps):
    '''
    finding the root by bitsection method
    get: polynom, start point of x, end point of x, epsilon  distance from root
    return: the root with distance of epsilon from the real root, number of iteration
    if no root was found return none 
    '''
    x = sympy.symbols('x')
    count = 0
    a = small_x
    b = large_x
    
    while ((b-a)/2 > eps):
        count = count + 1
        c = (a+b)/2
        if polynom.subs(x,c) == 0: # we found c as root
            return c, count
        elif (polynom.subs(x,a) * polynom.subs(x,c) < 0): #root is between a and c
            b = c
        elif (polynom.subs(x,c) * polynom.subs(x,b) < 0): #root is between c and b
            a = c 
        else: #no root found
            print("no root found")
            return None        
    #root is between a and b but it smaller that epsilon
    return (a+b) /2, count


def Newton_Raphson(polynom, small_x, large_x, eps):
    '''
    finding the root by secant method
    get: polynom, start point of x, end point of x, epsilon  distance from root
    return: the root with distance of epsilon from the real root, number of iteration
    if derivative in x0 = 0 divide by zero and root was found return none 
    '''
    x = sympy.symbols('x')
    count = 0
    derivative = sympy.diff(polynom, x)
    x0 = small_x    
    
    while True:
        count = count + 1
        y0 = polynom.subs(x,x0)
        y_tag = derivative.subs(x,x0)
        
        if (y_tag == 0):
            print("derivative was 0 cant divide by 0 no root found")
            return None
        x1 = x0 - (y0 / y_tag)
        if abs(x1 - x0) < eps:
            return x1, count
        x0 = x1
        if count > 1000:
            print("Eror got more than 1000 iteration") 
            return None
        

def secant(polynom, small_x, large_x, eps):
    '''
    finding the root by secant method
    get: polynom, start point of x, end point of x, epsilon  distance from root
    return: the root with distance of epsilon from the real root, number of iteration
    if y0 = y1 we get divide by zero and root was found return none 
    '''
    x = sympy.symbols('x')
    count = 0
    x0 = small_x
    x1 = large_x
    y0 = polynom.subs(x,x0)
    
    while True:
        count = count + 1        
        y1 = polynom.subs(x,x1)
        if y1-y0 == 0: #evoid zero divide
            print("divide by zero cant continue")
            return None
        x2 = x1 - y1 * (x1-x0)/(y1-y0)
        if abs(x2-x1) < eps: #we are in epsilon distance from root
            return x2, count
        x0 = x1
        y0 = y1
        x1 = x2
        if count > 1000:
            print("Eror got more than 1000 iteration") 
            return None
    

print("enter the method for root finding:")
print("1. for Bisection, 2. Newton-Raphson, 3. Secant")
method = int(input("enter number"))
epsilon = 0.0001
x = sympy.symbols('x')
polynom = x**4 - 4*x**3 +x**2 + 6*x #

if not (method in [1, 2, 3]):
    print("Wrong input, try again.")

minRange = -6.001
maxRange = 6
steps_dis = 0.1
num_of_steps = int((maxRange-minRange)/steps_dis)#steps_dis is float and num_of_steps need to be int for for loop
roots=[]
for i in range(num_of_steps):
    a = minRange + i *steps_dis
    b = minRange + (i+1) * steps_dis
    
    ya = polynom.subs(x,a)
    yb = polynom.subs(x,b)
    if ya == 0: #the steps limites are roots
        roots.append((a,0))
    if yb == 0 and b == maxRange: #skip all yb = 0 except the end point
        roots.append((b,0))
    if ya * yb < 0:
        if method == 1:
            roots.append(bitsection(polynom,a,b,epsilon))
        elif method == 2:
            roots.append(Newton_Raphson(polynom,a,b,epsilon))
        elif method == 3:
            roots.append(secant(polynom,a,b,epsilon))

print("roots that we found: ",roots)




