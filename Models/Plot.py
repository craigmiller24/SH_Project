from matplotlib.ticker import ScalarFormatter
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import requests

variant_names = {
    'Other' : 'Other',
    'V-20DEC-01 (Alpha)' : 'Alpha',
    'V-21APR-02 (Delta B.1.617.2)' : 'Delta B.1.617.2',
    'V-21OCT-01 (Delta AY 4.2)' : 'Delta AY 4.2',
    'V-22DEC-01 (Omicron CH.1.1)' : 'Omicron CH.1.1',
    'V-22JUL-01 (Omicron BA.2.75)' : 'Omicron BA.2.75',
    'V-22OCT-01 (Omicron BQ.1)' : 'Omicron BQ.1',
    'V-22OCT-02 (Omicron XBB)' : 'Omicron XBB',
    'V-23APR-01 (Omicron XBB.1.16)' : 'Omicron XBB.1.16',
    'V-23AUG-01 (Omicron BA.2.86)' : 'Omicron BA.2.86',
    'V-23JAN-01 (Omicron XBB.1.5)' : 'Omicron XBB.1.5',
    'V-23JUL-01 (Omicron EG.5.1)' : 'Omicron EG.5.1',
    'VOC-21NOV-01 (Omicron BA.1)' : 'Omicron BA.1',
    'VOC-22APR-03 (Omicron BA.4)' : 'Omicron BA.4',
    'VOC-22APR-04 (Omicron BA.5)' : 'Omicron BA.5',
    'VOC-22JAN-01 (Omicron BA.2)' : 'Omicron BA.2',
}
# Get Dataset 1 from online source
endpoint = (
    'https://api.coronavirus.data.gov.uk/v1/data?'
    'filters=areaType=nation;areaName=england&'
    'structure={"date":"date","newCases":"newCasesByPublishDate"}'
)

# Function to fetch data
def get_data(url):
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()
    else:
        raise ValueError(f"Failed to fetch data: {response.status_code}")

# Convert Data to a DataFrame
data = get_data(endpoint)['data']
df = pd.DataFrame(data)

# Format and sort data values
df['date'] = pd.to_datetime(df['date'], errors='coerce')
df.dropna(subset=['date', 'newCases'], inplace=True)
df.sort_values(by='date', inplace=True)

# Get dataset 2 from csv file, Convert to DataFrame and format
df2 = pd.read_csv('DataSets/Reduced_England_Variants.csv', names=['AreaCode', 'AreaName', 'AreaType', 'date', 'Variant', 'cumWeeklySequenced', 'newWeeklyPercentage'], header=None)
df2.drop(columns=['AreaCode', 'AreaName', 'AreaType', 'cumWeeklySequenced'], inplace=True)
df2['date'] = pd.to_datetime(df2['date'])
df2.sort_values(by='date', inplace=True)

# Plot set up
fig, ax = plt.subplots()
fig.set_figheight(8)
fig.set_figwidth(18)
ax.set_title("Coronavirus Variant Breakdown - England",fontweight='bold')
ax.set_xlabel("Time")

# Plot newCases vs Date on a secondary axis
ax.bar(df['date'], df['newCases'], color='grey', label='Daily Cases',width=2)
ax.set_ylabel("Daily Infections")

# Plot newWeeklyPercentage vs date for each variant
ax2 = ax.twinx()
for variant, group in df2.groupby('Variant'):
    if variant_names[variant] not in ['Alpha','Delta B.1.617.2','Omicron XBB.1.16','Omicron BA.1','Omicron EG.5.1','Omicron BA.4','Omicron XBB','Omicron BA.2','Other','Omicron BA.2.75']:
        ax2.plot(group['date'], group['newWeeklyPercentage'], label=variant_names[variant])

# Formatting
ax.set_ylim(0,190_000)
ax2.set_ylabel("Percentage of Cases by Variant")
ax2.set_ylim(0,100)
ax2.xaxis.set_major_formatter(mdates.DateFormatter('%d-%m-%Y'))
ax2.xaxis.set_tick_params(rotation=0)
ax2.set_xlim(min(df2['date']),max(df2['date']))

# Combine legends
lines, labels = ax.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax.legend(lines + lines2, labels + labels2, loc='upper center', bbox_to_anchor=(0.5, -0.075), fancybox=True, ncol=7)

plt.tight_layout()
plt.show()
