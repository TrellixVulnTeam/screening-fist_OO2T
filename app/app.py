import pandas as pd
import dash
from dash import Dash, html, dcc
import plotly.express as px

from dash.dependencies import Input, Output
from dash.exceptions import PreventUpdate


app = dash.Dash(
    __name__,
    meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}],
)

server = app.server

df = pd.read_csv('data/mol-props.csv')

fig = px.scatter(df, 
                 x="ExactMolWt", 
                 y="NumHBD",
                 hover_data=['Item Name', 'ExactMolWt', 'NumHBD'])

app.layout = html.Div(
    [
        html.Div([
            html.Img(src='assets/logo.png', 
                     className="app__banner",
                     height='128',
                     width='128'),
            html.H1('Screening Fist', 
                     className="title"),
    dcc.Graph(
        id='example-graph',
        figure=fig
    )
            ])])

if __name__ == "__main__":
    app.run_server(debug=True)
