import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import re
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
from matplotlib_venn import venn3, venn3_circles, venn3_unweighted
from bokeh.layouts import gridplot
from bokeh.models import BoxSelectTool, LassoSelectTool, ColumnDataSource, \
DataRange1d, Select, DataTable, ColumnDataSource, DataTable, DateFormatter, TableColumn
from bokeh.plotting import curdoc, figure
from bokeh.io import output_file, show
from bokeh.models import HoverTool
from bokeh.plotting import figure
from bokeh.models import Title
sns.set()


class move:
    '''A user-created class: `MOVE` object.
    Params:
    -----------

    Attributes:
    -----------
    intensity: dataframe
        The total intensity of each sample, subsets from the combined_protein.tsv 

    peptide: dataframe
        The cleaned peptide information from the combined_peptide.tsv file

    metadata: dataframe
        The summary information from the combined_protein.tsv, including gene names, 
        protein description, protein Existence etc. 

    samples: dataframe
        Container for sample type and replicates information.

    _log: dataframe
        The log2 transformed intensity.
        
    _median: dataframe
        The median normalized intensity.'''

    def __init__(self):

        # user accessible attributes
        self.intensity = pd.DataFrame()
        self.peptide = pd.DataFrame()
        self.metadata = pd.DataFrame()
        self.samples = pd.DataFrame()

        # attributes for manipulation and storage
        self._log = pd.DataFrame()
        self._median = pd.DataFrame()

    def _preprocess(self, protein_df, peptide_df, species, rename = None):
        '''Subsets a combined_protein dataframe protein_df to
        required columns and renames. Returns two parsed datafames:  
        Tabel 1 df_meta and Tabel 2 intensities df_int. 
        Exclude NaN Protein ID from the peptide_df. Returns the cleaned peptide dataframe.

        Params:
        -------
        protein_df: dataframe
            The protein dataframe to be prased

        peptide_df: dataframe
            The peptide dataframe to be cleaned

        species: string
            The name of the target species to subset the protein_df 

        rename: dict
            A dictionary to rename intensity columns in the protein_df. Defaults to
            None.

        Returns:
        --------
        df_meta: dataframe
            A parsed summary and identification table. 

        df_int: dataframe
            A parsed intensity table. It subsets to protein ID, total intensity,  
            and renames the columns.
        
        peptide_df: dataframe
            A cleaned peptide table. NaN proteind Ids excluded. '''

        # restrict organism to user input species 
        protein_df = protein_df[protein_df['Organism'].str.contains(species, flags=re.IGNORECASE, regex=True)]

        # set index as Protein ID
        protein_df.set_index('Protein ID',inplace=True)

        # columns needed for Tabel1 Protein summary and identification
        meta_cols= ['Protein Group',
                    'SubGroup',
                    'Protein',
                    'Entry Name',
                    'Gene Names',
                    'Protein Length',
                    'Coverage',
                    'Organism',
                    'Protein Existence',
                    'Description',
                    'Protein Probability',
                    'Top Peptide Probability',
                    'Unique Stripped Peptides',
                    'Summarized Total Spectral Count',
                    'Summarized Unique Spectral Count',
                    'Summarized Razor Spectral Count',
                    'Indistinguishable Proteins']

        # output Tabel 1 
        df_meta = protein_df[meta_cols]
        df_meta.to_csv("Table1_Protein_identification.csv")

        # create Table 2 Protein intensity
        df_int = protein_df[protein_df.columns[protein_df.columns.str.contains('Total Intensity')]]
        df_int = df_int.rename(columns = lambda x: x.replace('Total Intensity','').strip())

        # rename the intensity columns using user defined sample names
        if rename is not None:
            for col in df_int.columns:
                if col not in rename:
                    print(f"Can't rename {col}, plases fix your sample file")
            df_int = df_int.rename(columns = rename)

        # output Tabel 2
        df_int.to_csv("Table2_Protein_intensity.csv")

        # subset columns need from peptide dataframe
        peptide_cols = ['Sequence',
                        'Charge States',
                        'Probability',
                        'Assigned Modifications',
                        'Gene',
                        'Protein',
                        'Protein ID',
                        'Protein Description']
        peptide_df = peptide_df[peptide_cols]
        
        # For some sequence, the Protein ID is NaN and Probalibity is 0, need to be dropped
        peptide_df = peptide_df.dropna(subset=['Protein ID'])
        peptide_df.set_index('Protein ID',inplace=True)

        return df_meta, df_int, peptide_df

    def add_dataframes(self, protein_df, sample_df, peptide_df, species = 'Homo sapiens'):
        '''Add dataframes to the MOVE object.

        Params:
        -------
        protein_df: dataframe
            The dataframe of a combined_protein dataframe after pd.read_csv()

        protein_df: dataframe
            The dataframe of a peptide dataframe after pd.read_csv()

        sample_df: dataframe
            The dataframe of the user_sample_names.xlsx file after pd.read_excel()

        species: string
            The name of user input species to subset the protein_df to. Defaults to
            Homo sapiens.
            
        Returns:
        --------
        None'''

        self.samples = sample_df
        # create rename dict to store sample type and sample names
        rename = dict(zip(self.samples.iloc[:,0], self.samples.iloc[:,1]))
        # preprocess protein_df and peptide_df
        self.metadata, self.intensity, self.peptide = self._preprocess(protein_df, peptide_df, species = species, rename = rename)
        return None

    def log2_transform(self, figsize= (10,10)):
        '''Log2 tranform the intensity attribute of the MOVE object.

        Params:
        -------
        figsize (float, float), optional, default: (10,10)
            width, height of log transformed intensity distribution figure in inches. 
            
        Returns:
        --------
        None'''

        self._log = self.intensity.replace(0,np.nan)
        self._log = np.log2(self._log)

        # visualize 
        plt.figure(figsize=figsize)
        sns.boxplot(data=self._log, orient="h")

        # output figure and log transformed dataframe
        self._log.to_csv("D1_log2.csv")
        plt.savefig("D1_log.png", dpi=300, bbox_inches='tight')
        return None

    def median_normalize(self, figsize= (10,10)):
        '''Median normalize the _log intensity attribute of the MOVE object.

        Params:
        -------
        figsize (float, float), optional, default: (10,10)
            width, height of log transformed intensity distribution figure in inches. 
            
        Returns:
        --------
        None'''

        self._median = self._log - self._log.median()

        # visualize 
        plt.figure(figsize=figsize)
        sns.boxplot(data=self._median, orient="h")

        # output figure and log transformed dataframe
        self._median.to_csv("D2_median.csv")
        plt.savefig("D2_median.png", dpi=300, bbox_inches='tight')
        return None

    def _find_common_protein(self, threshold, sample_types):
        '''Select common proteins - subset of proteins with more than the threshold percent observed in at least one of the sample groups.

        Params:
        -------
        threshold: float
            The percent number to select common proteins
        
        sample_types: dict
            The mapping dict of sample type to each column 
            
        Returns:
        --------
        A set of common protein names'''

        protein_list = []
        for sample_type, col_names in sample_types.items():
            # calculate the percent observed for each protein
            percent = self._median[col_names].count(axis=1)/len(col_names)
            # select the common proteins (obversed more than the threshold)
            protein_list.extend(percent[percent >= threshold].index.tolist())
        return set(protein_list)
    
    def noise_filter(self, threshold = 0.65):
        '''Fill out noise from the _median intensity attribute of the MOVE object.

        Params:
        -------
        threshold, float, optional, default: 0.65
            The percent number to cut off noise proteins/select common proteins
        
        Returns:
        --------
        df_common: dataframe
            The dataframe of the selected common proteins'''

        # create the mapping dictionary of the sample type and column name
        sample_types = {}
        for sample_type in self.samples['Sample type'].unique():
            sample_types[sample_type] =  self.intensity.columns[self.intensity.columns.str.contains(sample_type)].tolist()
        
        # select common proteins and subset a common protein dataframe to output
        common_proteins = self._find_common_protein(threshold = threshold, sample_types = sample_types)
        df_common = self._median[self._median.index.isin(common_proteins)]
        df_common.to_csv("D3_common_proteins.csv")
        return df_common
    
    def plot_common_proteins(self, df_common, figsize = (15,8)):
        '''Create bar plot to visualize common proteins identified in each sample

        Params:
        -------
        df_common: dataframe
            The identified common proteins

        figsize (float, float), optional, default: (15,8)
            width, height of the bar plot figure in inches. 
        
        Returns:
        --------
        None'''

        pro_count = df_common.count()

        # assign bar color to each sample type
        samples = dict(zip(self.samples.iloc[:,1], self.samples.iloc[:,2]))
        cat_pal = sns.color_palette("pastel",self.samples['Sample type'].nunique())
        sample_colors = dict(zip(self.samples['Sample type'].unique(), cat_pal))
        colors = []
        for sample in pro_count.index:
            colors.append(sample_colors[samples[sample]])
        
        # plot the figure
        plt.figure(figsize=figsize)
        ax = pro_count.plot(kind='bar',color=colors)
        ax.set_ylabel('Protein Count',size ='large')
        plt.title("Identified common proteins in each sample",fontsize='large')
        plt.savefig("Common_protein.png", dpi=300, bbox_inches='tight')

        return None

    def _find_identified(self, df_common, sample_type):
        '''Subset identified common proteins in specified sample type

        Params:
        -------
        df_common: dataframe
            The identified common proteins
            
        sample_type: string
            The specified sample type name

        Returns:
        --------
        A set of identified common protein in the specified sample type'''

        df_sample = df_common[df_common.columns[df_common.columns.str.contains(sample_type)]]
        count = df_sample.count(axis = 1)
        return set(df_sample[count > 0].index)

    def _check_prevalence(self, df_common, proteins, sample_type):
        '''Subset a set of unique proteins in specified sample type, and calculate its average expression, count, and prevalence

        Params:
        -------
        df_common: dataframe
            The identified common proteins

        proteins: set
            The set of unique proteins
            
        sample_type: string
            The specified sample type name

        Returns:
        --------
        df_sample: dataframe
            A dataframe containing the unique proteins prevalence information in the specified sample type'''
            
        # subset to the unique proteins and the sample type
        samples = df_common.columns[df_common.columns.str.contains(sample_type)]
        df_sample = df_common[samples]
        df_sample = df_sample[df_sample.index.isin(proteins)]

        # calculate prevalence information
        df_sample['AveExpr'] = df_sample[samples].mean(axis = 1)
        df_sample['Count'] = df_sample[samples].count(axis = 1)
        df_sample['%Prevalence'] = round(df_sample['Count']/len(samples)*100,2)

        return df_sample

    def unique_proteins(self, df_common, control, case):
        '''Subset identified common proteins in specified sample type

        Params:
        -------
        df_common: dataframe
            The identified common proteins
            
        control: string
            The sample type name of the control group

        case: string
            The sample type name of the case group

        Returns:
        --------
        df_unique: dataframe
            A dataframe containing the unique proteins prevalence, summary and identification information in the specified sample type'''

        control_set = self._find_identified(df_common = df_common, sample_type = control)
        case_set = self._find_identified(df_common = df_common, sample_type = case)

        # visualize by Venn Diagram
        #plt.figure(figsize=(4,4))
        v = venn2(subsets = (control_set,case_set), set_labels = (control, case),
            set_colors=('salmon', 'skyblue'), alpha = 0.8)
        v.get_label_by_id('A').set_x(-0.25)
        v.get_label_by_id('B').set_x(0.25)
        for text in v.set_labels:
            text.set_fontsize(13)
        plt.annotate('Uniquely expressed proteins', xy=v.get_label_by_id('010').get_position() + np.array([0, 0.05]), xytext=(70,70),
             ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
             arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='gray'))
        plt.savefig(f"venn_{control}_{case}.png", dpi=300, bbox_inches='tight')


        # create case group unique proteins
        unique_in_case = case_set.difference(control_set)
        df_unique = self._check_prevalence(df_common=df_common, proteins=unique_in_case, sample_type=case)

        # merge with summary identification data 
        df_unique = df_unique.merge(self.metadata, left_index=True, right_index=True)
        df_unique.to_csv(f"uniquely_expressed_in_{case}.csv")
        return df_unique

    @classmethod
    def _draw_dashboard(cls, source, case):
        '''Create Bokeh html dashboard to visualize the uniquely expressed proteins

        Params:
        -------
        source: dataframe
            The dataframe of the uniquely expressed proteins in the case group
            
        case: string
            The sample type name of the case group

        Returns:
        --------
        None'''

        TOOLS="pan,wheel_zoom,box_select,lasso_select,reset"

        # create the scatter plot
        p = figure(tools=TOOLS, width=600, height=600, min_border=10, min_border_left=50,
           toolbar_location="above", 
           x_axis_location=None, 
           y_axis_location=None,
           title="Peptide identification",
          #title_location="left"
          )
        p.background_fill_color = "#fafafa"
        p.select(BoxSelectTool).select_every_mousemove = False
        p.select(LassoSelectTool).select_every_mousemove = False

        # configure visual properties on a plot's title attribute
        #p.title.text_color = "orange"
        p.title.text_font_size = "18px"
        from bokeh.models import Title
        # add extra titles with add_layout(...)
        n = len(source)
        t = f'{n} uniquely expressed proteins in {case}'
        p.add_layout(Title(text=t, 
                    align="center",
                    text_color = "orange",
                    text_font_size = "15px"), "below")

        r = p.scatter('Peptide_count', 'Probability', size='size', 
                    #size = 10,
                    #color="#3A5785", 
                    color = 'Color',
                    alpha=0.6,
                    hover_color="lightgreen", hover_alpha=0.8,
                    source=source)

        p.add_tools(HoverTool(
        tooltips=[#("index", "$index"),
                ("('Peptide_count','Probability')", "($x, $y)"),
                ("Protein", "@ProID"),
                ("Number of peptides","@Peptide_count"),
                ("Gene", "@Gene")
                ],
        mode="mouse", point_policy="follow_mouse", renderers=[r]
        ))

        # create the horizontal histogram
        hhist, hedges = np.histogram(source['Peptide_count'], bins=20)
        hzeros = np.zeros(len(hedges)-1)
        hmax = max(hhist)*1.1

        LINE_ARGS = dict(color="#3A5785", line_color=None)

        ph = figure(toolbar_location=None, width=p.width, height=200, x_range=p.x_range,
                #y_range=(-hmax, hmax), 
                title = 'Peptide counts', #Number of peptide
                min_border=10, min_border_left=50, y_axis_location="right")
        ph.xgrid.grid_line_color = None
        ph.yaxis.major_label_orientation = np.pi/4
        ph.background_fill_color = "#fafafa"

        ph.quad(bottom=0, left=hedges[:-1], right=hedges[1:], top=hhist, color="white", line_color="#3A5785")
        hh1 = ph.quad(bottom=0, left=hedges[:-1], right=hedges[1:], top=hzeros, alpha=0.5, **LINE_ARGS)
        hh2 = ph.quad(bottom=0, left=hedges[:-1], right=hedges[1:], top=hzeros, alpha=0.1, **LINE_ARGS)

        # create the vertical histogram
        vhist, vedges = np.histogram(source['Probability'], bins=20)
        vzeros = np.zeros(len(vedges)-1)
        vmax = max(vhist)*1.1

        pv = figure(toolbar_location=None, width=200, height=p.height, 
                #x_range=(-vmax, vmax),
                title = 'Probability (mean)',
                y_range=p.y_range, min_border=10, y_axis_location="right")
        pv.ygrid.grid_line_color = None
        pv.xaxis.major_label_orientation = np.pi/4
        pv.background_fill_color = "#fafafa"

        pv.quad(left=0, bottom=vedges[:-1], top=vedges[1:], right=vhist, color="white", line_color="#3A5785")
        vh1 = pv.quad(left=0, bottom=vedges[:-1], top=vedges[1:], right=vzeros, alpha=0.5, **LINE_ARGS)
        vh2 = pv.quad(left=0, bottom=vedges[:-1], top=vedges[1:], right=vzeros, alpha=0.1, **LINE_ARGS)

        # adding the table information
        columns = [
                    TableColumn(field="ProID", title="Protein ID"),
                    TableColumn(field="Gene", title="Gene"),
                    TableColumn(field="Peptide_count", title="Number of peptide"),
                    TableColumn(field="Protein_prevalence", title="%Prevalence"),
                    TableColumn(field="%Probability", title="%Probability (mean)"),
                ]
        data_table = DataTable(source=ColumnDataSource(source), columns=columns, width=800, height=600)
        layout = gridplot([[p, pv, data_table], [ph, None]], merge_tools=False)

        # save
        output_file(f"{case}_peptide_verification.html")
        show(layout)
        
        return None

    def peptide_verification(self, case, df_unique, threshold = 50):
        '''Query the peptide level idenfication information for the uniquely expressed proteins in the user specified case group

        Params:
        -------
        df_unique: dataframe
            The uniquely expressed proteins in the case group
            
        case: string
            The sample type name of the case group

        threshold, float, default: 50
            The proteins which prevalence greater than or equal to the threshold will be highlighted as pink in the dashboard


        Returns:
        --------
        None'''

        # subsets a peptide dataframe for the uniquely expressed proteins
        df_p = self.peptide[self.peptide.index.isin(df_unique.index)]
        
        # create a summary dataframe, and calculate the number of unique peptide and the mean of the probability of those peptide
        df_sum = pd.concat([df_p.groupby(df_p.index).Sequence.nunique(), df_p.groupby(df_p.index).Probability.mean()],axis=1)
        df_sum.rename(columns = {'Sequence': 'Peptide_count'},inplace=True)

        # merge with protein's data
        df_sum = df_sum.merge(df_unique[['Gene Names','Protein','Description','Count','%Prevalence','Indistinguishable Proteins']], 
              left_index= True, right_index=True)
        df_sum = df_sum.reset_index()

        # rename for later locate the columns in the hovering tool
        df_sum.rename(columns= {'Count': 'Protein_count', '%Prevalence': 'Protein_prevalence','Protein ID': 'ProID',
                        'Gene Names':'Gene'},inplace=True)
        
        # sort by prevalence
        df_sum = df_sum.sort_values(by = ['Protein_prevalence','Peptide_count'],ascending=False)
        
        # Targets the proteins with higher than threshold prevalence that the user might be intrested to see
        targets = df_sum[df_sum['Protein_prevalence'] >= threshold]['Gene'].to_list()
        df_sum['Color'] = df_sum.Gene.apply(lambda x: '#bf9593' if x in targets else '#3A5785')
        
        # round probability
        df_sum['%Probability'] = round(df_sum['Probability']*100,3)

        # define the size of circle in the dashboard 
        df_sum['size'] = df_sum['Protein_prevalence']**0.8
        self._draw_dashboard(df_sum,case)

        return None
    
        




