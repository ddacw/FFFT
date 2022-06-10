# https://www.kaggle.com/datasets/muthuj7/weather-dataset

import pandas as pd

filename = "weatherHistory.csv"
fileout = "temp.txt"

with open(fileout, 'w') as out, open(filename, 'r') as file:
    df = pd.read_csv(file)
    column = df['Temperature (C)']
    print(column.head)
    out.write(column.to_string(header=False, index=False))
