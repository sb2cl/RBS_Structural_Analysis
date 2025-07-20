from flask import Flask, render_template, request, jsonify
import json
import pandas as pd

from classes.RBS import RBS

app = Flask(__name__)

# Load the CSV file with the cluster data
df_clusters = pd.read_csv('data/clusters_with_stats.csv')
# Convert the SEQ_unique column from a string to a list of strings. Values are separated by a comma
df_clusters['SEQ_unique'] = df_clusters['SEQ_unique'].apply(lambda x: sorted(x.split(',')))

# Load the CSV file with sequences, levels and clusters
df_sequences = pd.read_csv('data/sequences_levels_clusters.csv')
# Add the stats of clusters to the df_sequences dataframe, merging on the CLUSTER column
df_sequences = df_sequences.merge(df_clusters[['CLUSTER', 'CORE REL EXPR_mean', 'CORE REL EXPR_std']],
    on='CLUSTER')
# Replace spaces in the colum names with underscores
df_sequences.columns = df_sequences.columns.str.replace(' ', '_')
# Round the values of expression mean and std
df_sequences['CORE_REL_EXPR_mean'] = df_sequences['CORE_REL_EXPR_mean'].round(4)
df_sequences['CORE_REL_EXPR_std'] = df_sequences['CORE_REL_EXPR_std'].round(4)

# Initialize the RBS class
rbs = RBS()

# Load the SD Core data with expression levels and statistics
with open('data/sd_core_data.json', 'r', encoding='utf8') as json_file:
    SEQUENCE_DATA = json.load(json_file)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/browse_sequences', methods=['GET', 'POST'])
def browse_sequences():
    context = {}

    return render_template('browse_sequences.html', context=context, app_data=SEQUENCE_DATA)

@app.route('/calculate', methods=['POST'])
def calculate():
    """
    Return the statistics for a given input string.

    This function is a Flask route that handles POST requests to the '/calculate' endpoint.
    It expects a JSON payload with an 'inputString' key.
    The function returns a JSON response with the statistics for the input string.

    Returns:
        A JSON response with the statistical values.
    """

    # Extract the input string from the request payload
    data = request.json
    input_string = data.get('inputString', '')

    seq = SEQUENCE_DATA.get(input_string, None)

    return jsonify(seq)

@app.route('/sequence_to_etr', methods=['GET', 'POST'])
def sequence_to_tir():
    context = {}
    results = None

    if request.method == 'POST':
        # Get the sequence from the form
        context['full_sequence'] = request.form['sequence'].upper()

        # Determine the sequence parts
        context['upstream'], context['core'], context['spacing'], context['expr'] = rbs.predict_spacing(
            context['full_sequence'])

        # Check if core sequence is present in the cores that belong to a cluster
        results = df_clusters[
            df_clusters['SEQ_unique'].apply(lambda x: context['core'] in x if isinstance(x, list) else False)]

    return render_template('sequence_to_etr.html', context=context, results=results)


@app.route('/etr_to_sequence', methods=['GET', 'POST'])
def etr_to_sequence():
    etr_value = None
    results = None

    if request.method == 'POST':

        try:
            etr_value = float(request.form['etr_value'])

            # Filter rows based on the range CORE REL EXPR mean ± CORE REL EXPR std
            # Note: Using copy to avoid SettingWithCopyWarning since we are adding a column later
            results = df_clusters[(df_clusters['CORE REL EXPR_mean'] - df_clusters['CORE REL EXPR_std'] <= etr_value) &
                                  (df_clusters['CORE REL EXPR_mean'] + df_clusters[
                                      'CORE REL EXPR_std'] >= etr_value)].copy()

            # Calculate the normalized distance using the mean and standard deviation
            results['distance_normalized'] = abs(results['CORE REL EXPR_mean'] - etr_value) / results[
                'CORE REL EXPR_std']

            # Sort by normalized distance, so the first result is the best
            results = results.sort_values(by='distance_normalized')

        except ValueError:

            results = None

    return render_template('etr_to_sequence.html', etr_value=etr_value, results=results)


@app.route('/browse_clusters', methods=['GET'])
def browse_clusters():

    return render_template('browse_clusters.html', data=df_sequences)


if __name__ == '__main__':
    app.run(debug=True)
