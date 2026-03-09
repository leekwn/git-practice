import pandas as pd

# Make a dictionary (columns)
data = {
    "Name": ["Alice", "Bob", "Charlie"],
    "Age": [25, 32, 18],
    "City": ["Seoul", "Busan", "Incheon"]
}

# Convert to DataFrame
df = pd.DataFrame(data)

print(df)
