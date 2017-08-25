


def l_func():
    lst = []
    for i in range(0,3):
        lst.append(i*i)
    return lst

def y_func():
    print("inside yield")
    for i in range(0,3):
        yield i*i

gen = y_func()
gen.__next__()
for i in gen:
    print(i)

# It has been reaped/harvested!
for i in gen:
    print(i)

lst = l_func()

for i in lst:
    print (i)

for i in lst:
    print (i)