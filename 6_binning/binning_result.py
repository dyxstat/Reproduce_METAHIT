import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Data for bar charts
methods = ['METAHIT', 'bin3C', 'MetaCC', 'ImputeCC']

# Group 5 data (blue color scheme)
greater_90_human_gut_5 = np.array([69, 44, 49, 59])[::-1]  
greater_70_human_gut_5 = np.array([107, 67, 73, 85])[::-1]  
greater_50_human_gut_5 = np.array([121, 75, 82, 97])[::-1]  

greater_90_pig_gut_5 = np.array([42, 13, 34, 38])[::-1]  
greater_70_pig_gut_5 = np.array([86, 29, 57, 71])[::-1] 
greater_50_pig_gut_5 = np.array([108, 37, 68, 92])[::-1]

greater_90_bovine_skin_5 = np.array([37, 20, 34, 35])[::-1]  
greater_70_bovine_skin_5 = np.array([61, 31, 50, 51])[::-1]  
greater_50_bovine_skin_5 = np.array([74, 38, 59, 64])[::-1]  

greater_90_wastewater_5 = np.array([62, 20, 55, 57])[::-1]  
greater_70_wastewater_5 = np.array([151, 52, 105, 124])[::-1]  
greater_50_wastewater_5 = np.array([198, 76, 135, 165])[::-1] 

greater_90_mats_5 = np.array([13, 5, 8, 11])[::-1] 
greater_70_mats_5 = np.array([48, 16, 21, 35])[::-1] 
greater_50_mats_5 = np.array([88, 29, 38, 69])[::-1] 

greater_90_sheep_gut_5 = np.array([487, 250, 256, 390])[::-1]  
greater_70_sheep_gut_5 = np.array([698, 324, 343, 575])[::-1]  
greater_50_sheep_gut_5 = np.array([834, 365, 394, 700])[::-1]

# Group 10 data (red color scheme)
greater_90_human_gut_10 = np.array([70, 49, 65, 68])[::-1]  
greater_70_human_gut_10 = np.array([112, 72, 96, 96])[::-1]  
greater_50_human_gut_10 = np.array([131, 80, 108, 110])[::-1]  

greater_90_pig_gut_10 = np.array([46, 13, 57, 44])[::-1]  
greater_70_pig_gut_10 = np.array([94, 31, 88, 80])[::-1]  
greater_50_pig_gut_10 = np.array([117, 40, 101, 104])[::-1]  

greater_90_bovine_skin_10 = np.array([41, 21, 42, 44])[::-1] 
greater_70_bovine_skin_10 = np.array([67, 32, 61, 63])[::-1]  
greater_50_bovine_skin_10 = np.array([80, 39, 70, 76])[::-1]  

greater_90_wastewater_10 = np.array([81, 28, 78, 79])[::-1]  
greater_70_wastewater_10 = np.array([185, 65, 146, 162])[::-1]  
greater_50_wastewater_10 = np.array([238, 91, 182, 209])[::-1]  

greater_90_mats_10 = np.array([21, 6, 19, 24])[::-1] 
greater_70_mats_10 = np.array([74, 24, 45, 61])[::-1]  
greater_50_mats_10 = np.array([124, 39, 68, 99])[::-1]

greater_90_sheep_gut_10 = np.array([540, 282, 356, 499])[::-1]  
greater_70_sheep_gut_10 = np.array([771, 360, 466, 736])[::-1] 
greater_50_sheep_gut_10 = np.array([929, 403, 532, 881])[::-1]  

# Create subplots - 6 rows, 2 columns
fig = make_subplots(
    rows=6, cols=2,
    vertical_spacing=0.07,
    horizontal_spacing=0.14
)

# Environment labels for y-axis
environments = ["Human Gut", "Pig Gut", "Bovine Skin", "Wastewater", "Mats", "Sheep Gut"]

# Color schemes
blue_colors = ['#4F74B7', '#8198C0', '#C7D2E3']
red_colors = ['#DD514A', '#E3918F', '#F0C2C1']

# Group 5 data (left column)
group5_data = [
    (greater_90_human_gut_5, greater_70_human_gut_5, greater_50_human_gut_5),
    (greater_90_pig_gut_5, greater_70_pig_gut_5, greater_50_pig_gut_5),
    (greater_90_bovine_skin_5, greater_70_bovine_skin_5, greater_50_bovine_skin_5),
    (greater_90_wastewater_5, greater_70_wastewater_5, greater_50_wastewater_5),
    (greater_90_mats_5, greater_70_mats_5, greater_50_mats_5),
    (greater_90_sheep_gut_5, greater_70_sheep_gut_5, greater_50_sheep_gut_5)
]

# Group 10 data (right column)
group10_data = [
    (greater_90_human_gut_10, greater_70_human_gut_10, greater_50_human_gut_10),
    (greater_90_pig_gut_10, greater_70_pig_gut_10, greater_50_pig_gut_10),
    (greater_90_bovine_skin_10, greater_70_bovine_skin_10, greater_50_bovine_skin_10),
    (greater_90_wastewater_10, greater_70_wastewater_10, greater_50_wastewater_10),
    (greater_90_mats_10, greater_70_mats_10, greater_50_mats_10),
    (greater_90_sheep_gut_10, greater_70_sheep_gut_10, greater_50_sheep_gut_10)
]

# Add Group 5 plots (left column - blue)
for i, (data_90, data_70, data_50) in enumerate(group5_data):
    row = i + 1
    show_legend = (i == 0)
    
    fig.add_trace(
        go.Bar(y=methods[::-1], x=data_50, orientation='h', name="Comp ≥ 50%", 
               marker_color=blue_colors[2], showlegend=show_legend, 
               legendgroup="group1", legendgrouptitle_text="Cont < 5%"),
        row=row, col=1
    )
    fig.add_trace(
        go.Bar(y=methods[::-1], x=data_70, orientation='h', name="Comp ≥ 70%", 
               marker_color=blue_colors[1], showlegend=show_legend, 
               legendgroup="group1"),
        row=row, col=1
    )
    fig.add_trace(
        go.Bar(y=methods[::-1], x=data_90, orientation='h', name="Comp ≥ 90%", 
               marker_color=blue_colors[0], showlegend=show_legend, 
               legendgroup="group1"),
        row=row, col=1
    )

# Add Group 10 plots (right column - red)
for i, (data_90, data_70, data_50) in enumerate(group10_data):
    row = i + 1
    show_legend = (i == 0)
    
    fig.add_trace(
        go.Bar(y=methods[::-1], x=data_50, orientation='h', name="Comp ≥ 50%", 
               marker_color=red_colors[2], showlegend=show_legend, 
               legendgroup="group2", legendgrouptitle_text="Cont < 10%"),
        row=row, col=2
    )
    fig.add_trace(
        go.Bar(y=methods[::-1], x=data_70, orientation='h', name="Comp ≥ 70%", 
               marker_color=red_colors[1], showlegend=show_legend, 
               legendgroup="group2"),
        row=row, col=2
    )
    fig.add_trace(
        go.Bar(y=methods[::-1], x=data_90, orientation='h', name="Comp ≥ 90%", 
               marker_color=red_colors[0], showlegend=show_legend, 
               legendgroup="group2"),
        row=row, col=2
    )

fig.update_layout(
    font=dict(family="Arial"),
    height=2000, width=1400,
    barmode='overlay', 
    margin=dict(l=200, r=200, t=50, b=100),
    template="simple_white",
    legend=dict(
        orientation="h",
        font=dict(size=25, family="Arial"),
        x=0.5,
        y=-0.05,
        xanchor="center",
        yanchor="top",
        grouptitlefont=dict(size=30, family="Arial")
    )
)

# Update axes for all subplots
for row in range(1, 7):
    for col in [1, 2]:
        fig.update_xaxes(title=dict(text="Number of bins", font=dict(size=25, family="Arial")), row=row, col=col)
        fig.update_xaxes(tickfont=dict(size=25, family="Arial"), row=row, col=col)
        fig.update_yaxes(tickfont=dict(size=25, family="Arial"), row=row, col=col)

# Add environment labels (A, B, C, D, E, F) on the left
letters = ["A", "B", "C", "D", "E", "F"]
for i, letter in enumerate(letters):
    if i == 0:
        yref = "y domain"
    else:
        yref = f"y{2*i+1} domain"
    
    fig.add_annotation(
        font=dict(size=45, family="Arial"),
        x=-0.15, y=0.95,
        xref="paper", yref=yref,
        text=letter,
        showarrow=False,
        textangle=0,  # Horizontal text
        xanchor="right",
        yanchor="top"
    )

fig.write_image("binning_result.pdf")
