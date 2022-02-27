import pandas as pd
import dash
from dash import Dash, html, dcc
import plotly.express as px

from dash.dependencies import Input, Output
from dash.exceptions import PreventUpdate

from utils import make_dash_table
app = dash.Dash(
    __name__,
    meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}],
)

server = app.server

df_ = pd.read_csv('data/mol-props-3d.csv')
cols = ['Item Name', 
        'CrippenDescriptors', 'LabuteASA', 'TPSA', 'ExactMolWt',
        'NumLipinskiHBD', 'NumLipinskiHBA', 'NumHBD', 'NumHBA',
        'NumRotatableBonds', 'NumRings', 'NumAromaticRings',
        'NumSaturatedRings', 'NumHeterocycles', 'NumAromaticHeterocycles',
        'NumAromaticCarbocycles', 'NumSaturatedHeterocycles',
        'NumSaturatedCarbocycles', 'NumAliphaticRings',
        'NumAliphaticHeterocycles', 'NumAliphaticCarbocycles', 'NumHeteroatoms',
        'NumAmideBonds', 'FractionCSP3', 'Chi0v', 'Chi1v', 'Chi2v', 'Chi3v',
        'Chi4v', 'Chi0n', 'Chi1n', 'Chi2n', 'Chi3n', 'Chi4n', 'HallKierAlpha',
        'Kappa1', 'Kappa2', 'Kappa3', 'NumSpiroAtoms', 'NumBridgeheadAtoms',
        'x','y','z', 'img']

df = df_[cols]
IMG = df['img'][0]
DRUG = df['Item Name'][0]

app.layout = html.Div([
                html.Div([
                    html.Div([
                        html.Img(src='assets/logo.png', 
                                 className="app__banner",
                                 height='128',
                                 width='128'),
                        html.H1('Screening Fist', 
                                 className="title")],
                                 style={'textAlign':'center',
                                       'font-family':'noto sans'}),
                dcc.Dropdown(id='label',options=df.columns,value='NumLipinskiHBD'),
                dcc.Graph(id='graph',
                          hoverData={"points": [{"pointNumber": 0}]}),
                html.Div([html.Img(id="chem_img",
                                 className="chem__img",
                                 height='128',
                                 width='128',
                                 src=IMG),
                          html.P(DRUG, id='chem_name',
                              style={'font-family':'noto sans'})],
                         className="chem__img__container")
            ])])

@app.callback(Output("graph", "figure"), 
              Input("label", "value"))
def figure(label):
    fig = px.scatter_3d(df, 
                     x='x', 
                     y='y',
                     z='z',
                     color=label,
                     hover_data=['Item Name'],
                     template='plotly_dark')
    fig.update_layout(
    hoverlabel=dict(
        bgcolor="gray",
        font_color='black',
        font_size=16,
        font_family="Noto Sans"
    )
)

    return fig

@app.callback([Output('chem_img','src'),
               Output('chem_name','children')],
               Input('graph', 'hoverData'))
def hover(hoverData):
    point_number = hoverData["points"][0]["pointNumber"]
    row = df.iloc[point_number,:]
    img = 'assets/' + row['img']
    return img, row['Item Name']


if __name__ == "__main__":
    app.run_server(debug=True)
