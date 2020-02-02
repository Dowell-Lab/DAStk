import DAStk.process_atac as p
import DAStk.differential_md_score as d
import DAStk.barcode_plot as b
import DAStk.ma_plot as m
import DAStk.tf_result_explanations as e
import DAStk.tf_intersect as i

def process_atac():
    p.main()

def differential_md_score():
    d.main()
    
def barcode_plot():
    b.main() 
    
def ma_plot():
    m.main()     
    
def tf_result_explanations():
    e.main()     

def tf_intersect():
    i.main()     
