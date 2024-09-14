from flask import Flask, render_template, request
import pandas as pd

from classes.RBS import RBS

app = Flask(__name__)

# Load the CSV file with the cluster data
df = pd.read_csv('data/clusters_with_stats.csv')
# Convert the SEQ_unique column from a string to a list of strings. Values are separated by a comma
df['SEQ_unique'] = df['SEQ_unique'].apply(lambda x: sorted(x.split(',')))

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
        context['full_sequence'] = request.form['sequence']

        # Determine the sequence parts
        context['upstream'], context['core'], context['spacing'], context['expr'] = rbs.predict_spacing(context['full_sequence'])

        # Get the basic expression level
        context['basic_expr'] = rbs.SEQS[context['core']]

        # Check if core sequence is present in the cores that belong to a cluster
        results = df[df['SEQ_unique'].apply(lambda x: context['core'] in x if isinstance(x, list) else False)]

    return render_template('sequence_to_tir.html', context=context, results=results)

@app.route('/tir_to_sequence', methods=['GET', 'POST'])
def tir_to_sequence():
    result = None
    if request.method == 'POST':
        try:
            value = float(request.form['value'])
            # Filter rows based on the range CORE REL EXPR mean ± CORE REL EXPR std
            result = df[(df['CORE REL EXPR mean'] - df['CORE REL EXPR std'] <= value) &
                        (df['CORE REL EXPR mean'] + df['CORE REL EXPR std'] >= value)]
        except ValueError:
            result = None
    return render_template('tir_to_sequence.html', result=result)

if __name__ == '__main__':
    app.run(debug=True)

