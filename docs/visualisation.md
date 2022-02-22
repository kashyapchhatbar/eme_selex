---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Data visualisation

> Analysis performed on ZFC4 and no protein control samples from [E-MTAB-11484](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-11484/) 

## Enrichment of relatively GC-rich 5-mers

```{code-cell} ipython3
:tags: [remove-input]
import pandas as pd
import plotly.io as pio
import plotly.express as px
import plotly.offline as py

melt_fractions_mean = pd.read_csv("data/melt_fractions_mean.tsv", sep="\t")

fig = px.scatter(melt_fractions_mean[melt_fractions_mean["AT"].isin([1,2,3])], 
                 facet_col="AT", x="cycle_", y="value", color="protein", 
                 error_y="value_std", template='simple_white',
                 category_orders={"protein": ["None", "ZFC4"],
                                  "AT": [1,2,3]},
                 color_discrete_sequence=["grey", "dodgerblue"],
                 width=800, height=400, custom_data=["kmer"],
                 labels={
                     "cycle_": "SELEX Cycle",
                     "value": "abundance",
                     "AT": "No. of A/T in kmer"
                 })
fig.update_traces(marker=dict(size=8, line=dict(color='DarkSlateGrey', width=1)),
                  selector=dict(mode='markers'), 
                  error_y=dict(thickness=1, width=2, color='DarkSlateGrey'),
                  hovertemplate="<br>".join([                      
                      "Abundance: %{y:.3f}",
                      "kmer: %{customdata[0]}",
                  ]))
fig.for_each_xaxis(lambda xaxis: xaxis.update(dict(
                      tickmode = 'array',
                      tickvals = [0,2.5,4.5,6.5],
                      ticktext = ['0', '1', '3', '6']
                  )))
fig.update_layout(
    hoverlabel=dict(
        bgcolor="white",
        font_size=16,
    ),
    yaxis_range=[-0.05,0.65]
)
fig.show()
```

## Visualise the enrichment of relatively AT-rich 5-mers

```{code-cell} ipython3
:tags: [remove-input]
fig = px.scatter(melt_fractions_mean[melt_fractions_mean["AT"].isin([3,4,5])], 
                 facet_col="AT", x="cycle_", y="value", color="protein", 
                 error_y="value_std", template='simple_white',
                 category_orders={"protein": ["None", "ZFC4"],
                                  "AT": [3,4,5]},
                 color_discrete_sequence=["grey", "dodgerblue"],
                 width=800, height=400, custom_data=["kmer"],
                 labels={
                     "cycle_": "SELEX Cycle",
                     "value": "abundance",
                     "AT": "No. of A/T in kmer"
                 })
fig.update_traces(marker=dict(size=8, line=dict(color='DarkSlateGrey', width=1)),
                  selector=dict(mode='markers'), 
                  error_y=dict(thickness=1, width=2, color='DarkSlateGrey'),
                  hovertemplate="<br>".join([                      
                      "Abundance: %{y:.3f}",
                      "kmer: %{customdata[0]}",
                  ]))
fig.for_each_xaxis(lambda xaxis: xaxis.update(dict(
                      tickmode = 'array',
                      tickvals = [0,2.5,4.5,6.5],
                      ticktext = ['0', '1', '3', '6']
                  )))
fig.update_layout(
    hoverlabel=dict(
        bgcolor="white",
        font_size=16,
    ),
    yaxis_range=[-0.05,0.65]
)

fig.show()
```