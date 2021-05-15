import os
import sys
import time

num_of_astr = 500
last_time = 100

os.system("rm *.out")
os.system("rm *.exe")
os.system("rm *.txt")

f1 = open("num_of_astr_config.txt",'a')
f2 = open("last_time_config.txt","a")
f1.write(str(num_of_astr))
f2.write(str(last_time))
#f.write("\n")
#f.write(str(last_time))

f1.close()
f2.close()


#os.system("gcc IV.c -o IV.exe -lm")
#os.system("IV.exe")
os.system("gcc IV.c -o IV.out -lm")
os.system("./IV.out")

start_time = time.time()

#os.system("gcc ode.c -o ode.exe -lm")
#os.system("ode.exe")
os.system("gcc ode.c -o ode.out -lm")
os.system("./ode.out")

end_time = time.time()

#os.system("gcc conn.c -o conn.exe")
#os.system("conn.exe")
print("time for", num_of_astr, "astroids, past", last_time, "years is", end_time - start_time)
print('c file done')

# os.system("python exp.py")
os.system("python3 hist.py")
# os.system("python conn.py")





# print('all done')
