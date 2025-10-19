import numpy as np
import csv
from io import StringIO
import pandas as pd
import os
import scipy
import re
from math import sqrt

dirpath = os.path.dirname(os.path.abspath(__file__)) + "/"
print(dirpath)

##spectra = input("Enter the name of the csv file that contains MEP spectra: ")
#mode_input = input("Reflective mode or Transmissive mode? ")
#mode = mode_input.lower()

print("Note: MEP names are case-sensitive!")
poly1 = input("Enter name of MEP 1: ")
poly2 = input("Enter name of MEP 2: ")

df = pd.read_csv(poly1 + ".csv", delimiter=",")
df2 = pd.read_csv(poly2 + ".csv", delimiter=",")
theo_df = pd.DataFrame()

theo_df["Wavelength(nm)"] = df["Wavelength(nm)"]

i = 0.05
j = 1.0 - i

column_names = []

while i <= 0.951:
    
    absorbance_values = i * df[poly1] + j * df2[poly2]
    
    column_name = f"{poly1} {round(i*100)} : {poly2} {round(j*100)}"
    column_names.append(column_name)
    
    theo_df[column_name] = absorbance_values
    #print("i=" + str(i) + ", j=" + str(j))
   
    i += 0.050
    j = 1.0 - i

theo_df.to_csv("theo_2MEP.csv", index=False, sep=",")

print("    theo_2MEP.csv was successfully created.")

#cont = input("Shall we continue?  ")

#if cont.lower() == "no":
#    print("Exiting the program.")
#    exit()  
#elif cont.lower() == "yes":
#    print("Continuing program.")

theo_l = float(input("Target L*?  "))
theo_a = float(input("Target a*?  "))
theo_b = float(input("Target b*?  "))


file = str(dirpath) + "theo_2MEP.csv"

with open(file) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    my_list = list(csv_reader)
    full_file = []
    for i in range(1,2046):
        list_f = list(map(float, my_list[i]))
        full_file.append(list_f)
    file_names = my_list[0]
    file_names.pop(0)
csv_file.close()

def lambda_avg(full_file, column):
    counter = 1
    color_lambda = int(full_file[0][0])
    total_data = []
    absorbances = []
    i = 0
    while i < len(full_file):
        if int(full_file[i][0]) == color_lambda:
            absorbances.append(full_file[i][column])
            i += 1
        else:
            avg_abs = pow(10,-sum(absorbances)/len(absorbances))
            total_data.append([color_lambda,avg_abs])
            absorbances = []
            absorbances.append(full_file[i][column])
            color_lambda = int(full_file[i][0])
            i += 1
    return total_data

final_data = lambda_avg(full_file, 1)

######################################################
############## Import Photop CSV         #############
######################################################

file = str(dirpath) + "photopic_ref.csv"
with open(file) as ref:
    csv_reader = csv.reader(ref, delimiter=',')
    my_list = list(csv_reader)
    photop_ref_file = []
    for i in range(1,201):
        list_f = list(map(float, my_list[i]))
        photop_ref_file.append(list_f)
ref.close()

######################################################
##############Clean Data for Photop and Scotop Trans #
######################################################

def clean_data(final_data, frequency, initial_lambda, final_lambda):
    i = 0
    even_data = []
    while i < len(final_data):
        if final_data[i][0] % frequency == 0:
            if final_data[i][0] < final_lambda and final_data[i][0] > initial_lambda:
                even_data.append(final_data[i])
            i += 1
        else:
            i += 1
    return even_data

photopic_data = clean_data(final_data, 2, 379, 781)

######################################################
############## Photopic Transmission                 #
######################################################

def photopic_transmission(photopic_data, photop_ref_file):
    top = 0
    k = 0
    i = 0
    while i < len(photopic_data)-1:
        k += photop_ref_file[i][0]*photop_ref_file[i][2] + pow(10,-20)
        top += photopic_data[i][1]*photop_ref_file[i][0]*photop_ref_file[i][2]
        i += 1
    inverted_k = 1/k
    value = top*inverted_k*100
    return value

photop_trans = photopic_transmission(photopic_data, photop_ref_file)


######################################################
############## Scotopic Transmission                 #
######################################################

def scotopic_transmission(photopic_data, photop_ref_file):
    top = 0
    k = 0
    i = 0
    while i < len(photopic_data)-1:
        k += photop_ref_file[i][0]*photop_ref_file[i][2] + pow(10,-20)
        top += photopic_data[i][1]*photop_ref_file[i][1]*photop_ref_file[i][2]
        i += 1
    inverted_k = 1/k
    value = top*inverted_k*100
    return value
scotop_trans = scotopic_transmission(photopic_data, photop_ref_file)


color_coor_data = clean_data(final_data, 10, 379, 771)

######################################################
############## Color Coordinate Stuff                #
######################################################

file = str(dirpath) + "color_coor_ref.csv"
with open(file) as ref:
    csv_reader = csv.reader(ref, delimiter=',')
    my_list = list(csv_reader)
    color_coor_ref = []
    for i in range(1,41):
        list_f = list(map(float, my_list[i]))
        color_coor_ref.append(list_f)
ref.close()


def color_coordinates(color_coor_data, color_coor_ref):
    i = 0
    X, Y, Z = 0, 0, 0
    while i < len(color_coor_data)-1:
        X += color_coor_data[i][1]*color_coor_ref[i][0]
        Y += color_coor_data[i][1]*color_coor_ref[i][1]
        Z += color_coor_data[i][1]*color_coor_ref[i][2]
        i += 1
    Tc = Y/1000
    x = X/(X+Y+Z)
    y = Y/(X+Y+Z)
    z = Z/(X+Y+Z)
    return Tc,x,y,z


Tc, x, y, z = color_coordinates(color_coor_data, color_coor_ref)

def LAB(color_coor_data, color_coor_ref):
    xn = 95.047
    yn = 100
    zn = 108.883
    i = 0
    X, Y, Z = 0, 0, 0
    while i < len(color_coor_data)-1:
        X += (color_coor_data[i][1]*color_coor_ref[i][0])/1000
        Y += (color_coor_data[i][1]*color_coor_ref[i][1])/1000
        Z += (color_coor_data[i][1]*color_coor_ref[i][2])/1000
        i += 1
    Q_vec = [X/xn, Y/yn, Z/zn]
    f_Q_vec = []
    for Q in Q_vec:
        if Q > pow(6.0/29.0,3):
            f_Q = pow(Q,(1.0/3.0))
            f_Q_vec.append(f_Q)
        else:
            f_Q = 841/108.0*Q + 4/29.0
            f_Q_vec.append(f_Q)
    l_star = 116*f_Q_vec[1]-16.0
    a_star = 500*(f_Q_vec[0]-f_Q_vec[1])
    b_star = 200*(f_Q_vec[1]-f_Q_vec[2])
    c_starab = sqrt(pow(a_star,2.0)+pow(b_star,2.0))
    return l_star, a_star, b_star, c_starab



######################################################
############## Neutrality Data Organization          #
######################################################


def range_avg(low_c, high_c, index_of_set, data_set):
    i = index_of_set
    x = 0
    list = []
    while x <= (high_c - low_c):
        list.append(data_set[i][1])
        x += 1
        i += 1
    average = sum(list)/len(list)
    return i, average

######################################################
############## Neutrality Calculation                #
######################################################

def neutrality_avg(final_data,Tc):
    color_lambda = int(final_data[0][0])
    total_data = []
    transmit = []
    low_c = 430
    high_c = 490
    counter = 0
    i = 0
    while int(final_data[i][0]) < low_c:
        i += 1
    while counter < 9:
        counter += 1
        i, trans = range_avg(low_c, high_c, i, final_data)
        i = i-31
        transmit.append([high_c,trans])
        low_c = high_c-30
        high_c = low_c+60
    weights = [5, 10, 10, 10, 10, 10, 10, 5, 1]
    temp_sum = 0.0
    i = 0
    while i < len(weights):
        value = 100*abs(1-(transmit[i][1]/(Tc/100)))
        temp_sum += value*weights[i]
        i += 1
    return temp_sum/sum(weights)


neutrality = neutrality_avg(final_data, Tc)

##########################################
## delta E calculation from theoretical ##
##########################################

def del_E(l_star, a_star, b_star):
    delta_E = float(sqrt((theo_l - l_star)**2 + (theo_a - a_star)**2 + (theo_b - b_star)**2))
    return delta_E



header = ["file_name","Tp", "Ts", "Neutrality", "x", "y", "c*ab","l*", "a*", "b*", "delta_E"]
## Create array for file
final_full_data = []
final_full_data.append(header)


q = 1
while q < len(full_file[0]):
    current_row = []
    current_row.append(file_names[q-1])
    final_data = lambda_avg(full_file, q)
    ## Calculate Tp
    photopic_data = clean_data(final_data, 2, 379, 781)
    photop_trans = photopic_transmission(photopic_data, photop_ref_file)
    ## Calculate Ts
    scotop_trans = scotopic_transmission(photopic_data, photop_ref_file)
    ### On to Neutrality Cleaning ###
    color_coor_data = clean_data(final_data, 10, 379, 771)
    ## Calculate x and y
    Tc, x, y, z = color_coordinates(color_coor_data, color_coor_ref)
    ## Calulate Neutrality
    neutrality = neutrality_avg(final_data, Tc)
    ## Calculate c*ab, l*, a*, b*
    l_star, a_star, b_star, c_starab =  LAB(color_coor_data, color_coor_ref)
    ## Calculate delta_E
    delta_E = del_E(l_star, a_star, b_star)

    current_row.append(photop_trans)
    current_row.append(scotop_trans)
    current_row.append(neutrality)
    current_row.append(x)
    current_row.append(y)
    current_row.append(c_starab)
    current_row.append(l_star)
    current_row.append(a_star)
    current_row.append(b_star)
    current_row.append(delta_E)
    final_full_data.append(current_row)
    q += 1

with open("output.csv", "w", newline='') as f:
    writer = csv.writer(f)
    writer.writerows(final_full_data)

# Read the CSV file 
df = pd.read_csv("output.csv")

# Find the row with the lowest delta_E value
min_delta_E_row = df.loc[df["delta_E"].idxmin()]

min_delta_E = min_delta_E_row['delta_E']

# Get the corresponding file_name
file_name = min_delta_E_row["file_name"]

#extract percent values from file_name
match = re.search(r'\b(\d{2})\b.*\b(\d{2})\b', file_name)

if match:
    percent_poly1 = float(match.group(1))
    percent_poly2 = float(match.group(2))

print(f"    Best starting point for first 2 MEPs:  " + file_name)
print(f"    Delta_E (theo to target): {min_delta_E:.2f}")




cont2 = input("Do you want to blend a third MEP?  ")

if cont2.lower() == "no":
    print("    Exiting the program. Goodbye!")
    exit()  
elif cont2.lower() == "yes":
    print("    Great! Let's continue.")


###############################
#### repeat for third MEP #####
###############################

spectra = "theo_2MEP"
poly3 = input("Enter the name of the third MEP: ")

blend = file_name

df = pd.read_csv(spectra + ".csv", delimiter=",")
df2 = pd.read_csv(poly3 + ".csv", delimiter=",")

theo3poly_df = pd.DataFrame()

theo3poly_df["Wavelength(nm)"] = df["Wavelength(nm)"]

i = 0.05
j = 1.0 - i

column_names = []

while i <= 0.951:
    
    absorbance_values = i * df[blend] + j * df2[poly3]
    
    column_name = f"{poly1} {round(percent_poly1*i,2)} : {poly2} {round(percent_poly2*i,2)} : {poly3} {round(j*100,2)}"
    column_names.append(column_name)
    
    theo3poly_df[column_name] = absorbance_values
    #print("i=" + str(i) + ", j=" + str(j))
   
    i += 0.050
    j = 1.0 - i

theo3poly_df.to_csv("theo3MEP.csv", index=False, sep=",")

print("    theo3MEP.csv was successfully created.")





file = str(dirpath) + "theo3MEP.csv"

with open(file) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    my_list2 = list(csv_reader)
    full_file = []
    for i in range(1,2046):
        list_f = list(map(float, my_list2[i]))
        full_file.append(list_f)
    file_names = my_list2[0]
    file_names.pop(0)
csv_file.close()


final_data = lambda_avg(full_file, 1)


######################################################
##############Clean Data for Photop and Scotop Trans #
######################################################

photopic_data = clean_data(final_data, 2, 379, 781)

######################################################
############## Photopic Transmission                 #
######################################################


photop_trans = photopic_transmission(photopic_data, photop_ref_file)


######################################################
############## Scotopic Transmission                 #
######################################################


scotop_trans = scotopic_transmission(photopic_data, photop_ref_file)


color_coor_data = clean_data(final_data, 10, 379, 771)

######################################################
############## Color Coordinate Stuff                #
######################################################

file = str(dirpath) + "color_coor_ref.csv"
with open(file) as ref:
    csv_reader = csv.reader(ref, delimiter=',')
    my_list = list(csv_reader)
    color_coor_ref = []
    for i in range(1,41):
        list_f = list(map(float, my_list[i]))
        color_coor_ref.append(list_f)
ref.close()


Tc, x, y, z = color_coordinates(color_coor_data, color_coor_ref)


######################################################
############## Neutrality Calculation                #
######################################################


neutrality = neutrality_avg(final_data, Tc)

##########################################
## delta E calculation from theoretical ##
##########################################



header = ["file_name","Tp", "Ts", "Neutrality", "x", "y", "c*ab","l*", "a*", "b*", "delta_E"]
## Create array for file
final_full_data = []
final_full_data.append(header)


q = 1
while q < len(full_file[0]):
    current_row = []
    current_row.append(file_names[q-1])
    final_data = lambda_avg(full_file, q)
    ## Calculate Tp
    photopic_data = clean_data(final_data, 2, 379, 781)
    photop_trans = photopic_transmission(photopic_data, photop_ref_file)
    ## Calculate Ts
    scotop_trans = scotopic_transmission(photopic_data, photop_ref_file)
    ### On to Neutrality Cleaning ###
    color_coor_data = clean_data(final_data, 10, 379, 771)
    ## Calculate x and y
    Tc, x, y, z = color_coordinates(color_coor_data, color_coor_ref)
    ## Calulate Neutrality
    neutrality = neutrality_avg(final_data, Tc)
    ## Calculate c*ab, l*, a*, b*
    l_star, a_star, b_star, c_starab =  LAB(color_coor_data, color_coor_ref)
    ## Calculate delta_E
    delta_E = del_E(l_star, a_star, b_star)

    current_row.append(photop_trans)
    current_row.append(scotop_trans)
    current_row.append(neutrality)
    current_row.append(x)
    current_row.append(y)
    current_row.append(c_starab)
    current_row.append(l_star)
    current_row.append(a_star)
    current_row.append(b_star)
    current_row.append(delta_E)
    final_full_data.append(current_row)
    q += 1

with open("output3MEP.csv", "w", newline='') as f:
    writer = csv.writer(f)
    writer.writerows(final_full_data)

# Read the CSV file 
df = pd.read_csv("output3MEP.csv")

# Find the row with the lowest delta_E value
min_delta_E_row = df.loc[df["delta_E"].idxmin()]

min_delta_E = min_delta_E_row['delta_E']

# Get the corresponding file_name
file_name = min_delta_E_row["file_name"]

print(f"    Best starting point for 3 MEPs:  " + file_name)
print(f"    Delta_E (theo to target): {min_delta_E:.2f}")
print(" ")
print("    Great work! Goodbye! :D")
