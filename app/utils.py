import pandas as pd
import dash
from dash import Dash, html, dcc
import plotly.express as px

from dash.dependencies import Input, Output
from dash.exceptions import PreventUpdate

def make_dash_table(selection, df):
    """ Return a dash defintion of an HTML table from a Pandas dataframe. """

    df_subset = df.loc[df["NAME"].isin(selection)]
    table = []

    for index, row in df_subset.iterrows():
        rows = []
        rows.append(html.Td([row["NAME"]]))
        rows.append(html.Td([html.Img(src=row["IMG_URL"])]))
        rows.append(html.Td([row["FORM"]]))
        rows.append(
            html.Td([html.A(href=row["PAGE"], children="Datasheet", target="_blank")])
        )
        table.append(html.Tr(rows))

    return table
