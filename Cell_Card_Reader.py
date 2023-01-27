"""

Retrieving information about specific cell types from Cell Cards


"""

#import urllib.request
import requests
import json


test_url = 'https://cellcards.org/podocyte.php'

data = requests.get(test_url)
data.raise_for_status()
print(data)
print(data.status_code)
print(data.headers)
print(data.text)

print(type(data.text))

print(len(data.text.split('<p>')))










