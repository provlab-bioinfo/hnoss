import os
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from datetime import datetime
from freyja import get_abundance, calc_rel_growth_rates, get_color_scheme

def make_dashboard(agg_df, meta_df, thresh, title, introText,
                   outputFn, headerColor, bodyColor, scale_by_viral_load,
                   config, lineage_info, nboots, serial_interval, days,
                   grThresh):
    df_ab_lin, df_ab_sum, dates_to_keep = get_abundance(agg_df, meta_df,
                                                        thresh,
                                                        scale_by_viral_load,
                                                        config, lineage_info)

    calc_rel_growth_rates(df_ab_lin.copy(deep=True), nboots,
                          serial_interval, outputFn,
                          daysIncluded=days, grThresh=grThresh)

    fig = go.Figure()

    default_color_lin = {
        11: px.colors.qualitative.Vivid,
        24: px.colors.qualitative.Dark24
    }
    color_lin = get_color_scheme(df_ab_lin,
                                 default_color_lin,
                                 config.get('Lineages'))
    for j, col in enumerate(df_ab_lin.columns):

        fig.add_trace(go.Scatter(
            x=df_ab_lin.index, y=df_ab_lin[col],
            hoverinfo='x+y',
            name=col,
            mode='markers+lines',
            hovertemplate="%{y:.1f}%",
            line=dict(width=0.5, color=color_lin[col]),
            visible=False,
            stackgroup='one'
        ))

    default_color_sum = {
        11: px.colors.qualitative.Pastel,
        24: px.colors.qualitative.Light24
    }
    color_sum = get_color_scheme(df_ab_sum,
                                 default_color_sum,
                                 config.get('VOC'))
    for j, col in enumerate(df_ab_sum.columns):
        fig.add_trace(go.Scatter(
            x=df_ab_sum.index, y=df_ab_sum[col],
            hoverinfo='x+y',
            name=col,
            mode='markers+lines',
            hovertemplate="%{y:.1f}%",
            line=dict(width=0.5, color=color_sum[col]),
            visible=True,
            stackgroup='one',
        ))
    # if needed, drop dates with missing viral load metadata
    meta_df = meta_df.set_index('sample_collection_datetime')
    if len(dates_to_keep) < meta_df.shape[0]:
        meta_df = meta_df.loc[dates_to_keep]
        df_ab_sum = df_ab_sum.loc[dates_to_keep]
        df_ab_lin = df_ab_lin.loc[dates_to_keep]
    fig.add_trace(go.Scatter(
        x=meta_df.index, y=meta_df['viral_load'],
        hoverinfo='x+y',
        mode='markers+lines',
        hovertemplate="%{y:.1f} copies/L",
        line=dict(width=1.25, color='blue'),
        visible=False,
    ))

    # load scaled abundances
    df_ab_lin_s = df_ab_lin.multiply(meta_df.viral_load,
                                     axis=0) / 100.

    default_color_lin_s = {
        11: px.colors.qualitative.Vivid,
        24: px.colors.qualitative.Dark24
    }
    color_lin_s = get_color_scheme(df_ab_lin_s,
                                   default_color_lin_s,
                                   config.get('Lineages'))
    for j, col in enumerate(df_ab_lin_s.columns):

        fig.add_trace(go.Scatter(
            x=df_ab_lin_s.index, y=df_ab_lin_s[col],
            hoverinfo='x+y',
            name=col,
            mode='markers+lines',
            hovertemplate="%{y:.1f} copies/L",
            line=dict(width=0.5, color=color_lin_s[col]),
            visible=False,
            stackgroup='one'
        ))

    df_ab_sum_s = df_ab_sum.multiply(meta_df.viral_load,
                                     axis=0) / 100.

    default_color_sum_s = {
        11: px.colors.qualitative.Pastel,
        24: px.colors.qualitative.Light24
    }
    color_sum_s = get_color_scheme(df_ab_sum,
                                   default_color_sum_s,
                                   config.get('VOC'))
    for j, col in enumerate(df_ab_sum_s.columns):
        fig.add_trace(go.Scatter(
            x=df_ab_sum_s.index, y=df_ab_sum_s[col],
            hoverinfo='x+y',
            name=col,
            mode='markers+lines',
            hovertemplate="%{y:.1f} copies/L",
            line=dict(width=0.5, color=color_sum_s[col]),
            visible=False,
            stackgroup='one'
        ))
    fig.update_layout(updatemenus=[dict(type="buttons",
                      direction='right',
                      active=0,
                      bgcolor='lightskyblue',
                      x=0.85,
                      y=1.07,
                      buttons=list([
                                   dict(label="Variants",
                                        method="update",
                                        args=[{
                                             "visible":
                                             [False] * df_ab_lin.shape[1] +
                                             [True] * df_ab_sum.shape[1] +
                                             [False] +
                                             [False] * df_ab_lin_s.shape[1] +
                                             [False] * df_ab_sum_s.shape[1]},
                                            {"yaxis": {"title":
                                             'VOC Prevalence',
                                                       "ticksuffix": '%',
                                                       "range": [0, 100]}}]),
                                   dict(label="Lineages",
                                        method="update",
                                        args=[{
                                              "visible":
                                              [True] * df_ab_lin.shape[1] +
                                              [False] * df_ab_sum.shape[1] +
                                              [False] +
                                              [False] * df_ab_lin_s.shape[1] +
                                              [False] * df_ab_sum_s.shape[1]},
                                              {"yaxis": {"title":
                                               'Lineage Prevalence',
                                                         "ticksuffix": '%',
                                                         "range": [0, 100]}}]),
                                   dict(label="Viral Load",
                                        method="update",
                                        args=[{
                                              "visible":
                                              [False] * df_ab_lin.shape[1] +
                                              [False] * df_ab_sum.shape[1] +
                                              [True] +
                                              [False] * df_ab_lin_s.shape[1] +
                                              [False] * df_ab_sum_s.shape[1]},
                                              {"yaxis": {"title":
                                                         'Virus copies/L',
                                                         "range":
                                                         [0,
                                                          meta_df['viral_load']
                                                          .max() * 1.1]}}]),
                                   dict(
                                       label="Load Scaled Variants",
                                       method="update",
                                       args=[{
                                              "visible":
                                              [False] * df_ab_lin.shape[1] +
                                              [False] * df_ab_sum.shape[1] +
                                              [False] +
                                              [False] * df_ab_lin_s.shape[1] +
                                              [True] * df_ab_sum_s.shape[1]},
                                             {"yaxis": {"title":
                                              'Variant copies/L',
                                                        "range":
                                                        [0,
                                                         meta_df['viral_load']
                                                         .max() * 1.1]}}]),
                                   dict(
                                       label="Load Scaled Lineages",
                                       method="update",
                                       args=[{
                                              "visible":
                                              [False] * df_ab_lin.shape[1] +
                                              [False] * df_ab_sum.shape[1] +
                                              [False] +
                                              [True] * df_ab_lin_s.shape[1] +
                                              [False] * df_ab_sum_s.shape[1]},
                                             {"yaxis": {"title":
                                                        'Lineage copies/L',
                                                        "range":
                                                        [0,
                                                         meta_df['viral_load']
                                                         .max() * 1.1]}}])]))],
                      template="plotly_white",
                      hovermode="x unified",
                      xaxis=dict(hoverformat="%B %d, %Y"),
                      legend=dict(yanchor="top",
                                  y=0.99,
                                  xanchor="right",
                                  x=1.1,
                                  itemsizing='constant'))

    fig.update_layout(yaxis_title='Variant Prevalence',
                      yaxis_ticksuffix="%",
                      yaxis_range=[0, 100.],
                      width=1000,
                      height=600)
    fig.update_layout(margin=dict(b=0, t=10))

    fig.update_xaxes(dtick="6.048e+8", tickformat="%b %d", mirror=True,
                     showline=False, ticks="", showgrid=False)
    # generate plot as a div, to combine with additional html/css
    fig.write_html("div-plot.html", full_html=False, default_width='50%',
                   config={'displaylogo': False, 'displayModeBar': False})
    # Generate a web page with the plot by placing it in the template.
    locDir = os.path.abspath(os.path.join(os.path.realpath(__file__),
                             os.pardir))
    webpage = open(os.path.join(locDir,
                                'data/dashboard_template.html'),
                   "r").read()
    webpage = webpage.replace("{title}", title)
    webpage = webpage.replace("{introText}", introText)
    webpage = webpage.replace("{plot}", open("div-plot.html", "r").read())
    webpage = webpage.replace("{lastUpdated}",
                              str(datetime.now().strftime("%b-%d-%Y %H:%M")))
    webpage = webpage.replace("{headerColor}",
                              headerColor)
    webpage = webpage.replace("{bodyColor}",
                              bodyColor)
    webpage = webpage.replace("{table}",
                              pd.read_csv(
                                  outputFn.replace('.html',
                                                   '_rel_growth_rates.csv'))
                                .to_html(index=False)
                                .replace('dataframe',
                                         'table table-bordered table-hover' +
                                         ' table-striped table-light' +
                                         ' table-bordered'))

    with open(outputFn, 'w') as outfile:
        outfile.write(webpage)
    os.remove('div-plot.html')
    print("Dashboard html file saved to " + outputFn)