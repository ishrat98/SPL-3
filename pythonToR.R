library(reticulate)

os <- import("os")
os$listdir(".")
os$listdir


py_install("pandas")
import pandas

Wrapper 
def read_flights(file):
  flights = pandas.read_csv(file)
flights = flights[flights['dest'] == "ORD"]
flights = flights[['carrier', 'dep_delay', 'arr_delay']]
flights = flights.dropna()
return flights