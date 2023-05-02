import argparse
import json
import pandas as pd
from collections import defaultdict

def parse_trace_file(trace_file):
    """
    Parses a Nextflow trace file and returns a Pandas DataFrame with the start and end times of each task.
    """
    # Read the trace file as a CSV file using Pandas
    df = pd.read_csv(trace_file,sep="\t")
    print(df)
    # Convert the 'time' column to datetime format
    df['time'] = pd.to_datetime(df['submit'])

    # Create a new DataFrame with the event data
    events = pd.DataFrame(columns=['start', 'end', 'label'])

    # Iterate over the rows of the DataFrame and extract the relevant data
    for i, row in df.iterrows():
        if row['type'] == 'TASK':
            label = f"{row['name']} [{row['id']}]"
            events = events.append({'start_time': row['time'], 'end_time': row['time'] + pd.Timedelta(seconds=row['duration']), 'label': label}, ignore_index=True)
    print(events)
    return pd.DataFrame.from_dict(events, orient='index')

def generate_timeline_report(timeline_df, output_file):
    """
    Generates an HTML timeline report using the Google Charts Timeline API.
    """
    # Define the colors to use for each task
    colors = ['red', 'green', 'blue', 'orange', 'purple', 'yellow', 'brown', 'gray']

    # Define the HTML code for the timeline report
    html = '<html>\n<head>\n'
    html += '  <script type="text/javascript" src="https://www.gstatic.com/charts/loader.js"></script>\n'
    html += '  <script type="text/javascript">\n'
    html += '    google.charts.load("current", {"packages":["timeline"]});\n'
    html += '    google.charts.setOnLoadCallback(drawChart);\n'
    html += '    function drawChart() {\n'
    html += '      var container = document.getElementById("timeline");\n'
    html += '      var dataTable = new google.visualization.DataTable();\n'
    html += '      dataTable.addColumn({ type: "string", id: "Task" });\n'
    html += '      dataTable.addColumn({ type: "string", id: "Label" });\n'
    html += '      dataTable.addColumn({ type: "date", id: "Start" });\n'
    html += '      dataTable.addColumn({ type: "date", id: "End" });\n'
    for i, task_name in enumerate(timeline_df.index):
        start_time = timeline_df.loc[task_name, 'start_time']
        end_time = timeline_df.loc[task_name, 'end_time']
        label = task_name
        color = colors[i % len(colors)]
        html += f'      dataTable.addRows([\n'
        html += f'        ["{task_name}", "{label}", new Date("{start_time}"), new Date("{end_time}")],\n'
        html += f'      ]);\n'
    html += '      var options = {\n'
    html += '        height: 1000,\n'
    html += '        colors: colors,\n'
    html += '      };\n'
    html += '      var chart = new google.visualization.Timeline(container);\n'
    html += '      chart.draw(dataTable, options);\n'
    html += '    }\n'
    html += '  </script>\n'
    html += '</head>\n<body>\n'
    html += '  <div id="timeline"></div>\n'
    html += '</body>\n</html>'

    # Write the HTML code to the output file
    with open(output_file, 'w') as f:
        f.write(html)


if __name__ == '__main__':
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Generate a timeline report from a Nextflow trace file')
    parser.add_argument('trace_file', help='the Nextflow trace file to parse')
    parser.add_argument('output_file', help='the output HTML file for the timeline report')
    args = parser.parse_args()

    # Parse the Nextflow trace file and generate the timeline report
    timeline_df = parse_trace_file(args.trace_file)
    generate_timeline_report(timeline_df, args.output_file)