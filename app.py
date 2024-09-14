from flask import Flask, render_template, request
import pandas as pd

from classes.RBS import RBS

app = Flask(__name__)

# Load the CSV file with the cluster data
df_clusters = pd.read_csv('data/clusters_with_stats.csv')
# Convert the SEQ_unique column from a string to a list of strings. Values are separated by a comma
df_clusters['SEQ_unique'] = df_clusters['SEQ_unique'].apply(lambda x: sorted(x.split(',')))

# Initialize the RBS class
rbs = RBS()


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/sequence_to_tir', methods=['GET', 'POST'])
def sequence_to_tir():
    context = {}
    results = None

    if request.method == 'POST':
        # Get the sequence from the form
        context['full_sequence'] = request.form['sequence'].upper()

        # Determine the sequence parts
        context['upstream'], context['core'], context['spacing'], context['expr'] = rbs.predict_spacing(
            context['full_sequence'])

        # Get the basic expression level
        context['basic_expr'] = rbs.SEQS[context['core']]

        # Check if core sequence is present in the cores that belong to a cluster
        results = df_clusters[
            df_clusters['SEQ_unique'].apply(lambda x: context['core'] in x if isinstance(x, list) else False)]

    return render_template('sequence_to_tir.html', context=context, results=results)


@app.route('/tir_to_sequence', methods=['GET', 'POST'])
def tir_to_sequence():
    tir_value = None
    results = None

    if request.method == 'POST':

        try:
            tir_value = float(request.form['tir_value'])

            # Filter rows based on the range CORE REL EXPR mean ± CORE REL EXPR std
            # Note: Using copy to avoid SettingWithCopyWarning since we are adding a column later
            results = df_clusters[(df_clusters['CORE REL EXPR_mean'] - df_clusters['CORE REL EXPR_std'] <= tir_value) &
                                  (df_clusters['CORE REL EXPR_mean'] + df_clusters[
                                      'CORE REL EXPR_std'] >= tir_value)].copy()

            # Calculate the normalized distance using the mean and standard deviation
            results['distance_normalized'] = abs(results['CORE REL EXPR_mean'] - tir_value) / results[
                'CORE REL EXPR_std']

            # Sort by normalized distance, so the first result is the ebst
            results = results.sort_values(by='distance_normalized')

        except ValueError:

            results = None

    return render_template('tir_to_sequence.html', tir_value=tir_value, results=results)


if __name__ == '__main__':
    app.run(debug=True)
