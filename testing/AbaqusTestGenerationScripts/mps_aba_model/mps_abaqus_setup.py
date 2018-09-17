# Build a 1x1x1 cube for testing mp simulation versus Abaqus
from mesh           import *
from interaction    import *
from sketch         import *
from part           import *
from material       import *
from section        import *
from assembly       import *
from step           import *
from load           import *
from mesh           import *
from job            import *

import mps_abaqus_val as val
reload(val)

def setup_model():
    my_model = mdb.models.values()[0]

    # Create sketch 
    my_sketch = my_model.ConstrainedSketch(
        name='__profile__', sheetSize=20)
    my_sketch.rectangle(point1=(0.0, 0.0), point2=(1.0, 1.0))

    
    # Create the solid/shell
    my_part = my_model.Part(name='CUBE', dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    my_part.BaseSolidExtrude(sketch=my_sketch, depth=1.0)
    del my_sketch

    # Create the material
    my_material = create_material(my_model)

    # Section
    my_section = my_model.HomogeneousSolidSection(name='Section-1', 
        material=val.MTRL_NAME, thickness=None)
    
    my_cell = my_part.cells.findAt((( 0.5, 0.5, 0.5),),)
    my_part.SectionAssignment(
        region=Region(cells = my_cell), 
        sectionName='Section-1', 
        thicknessAssignment=FROM_SECTION)
                
    # Create mesh
    create_mesh(my_part, my_cell)

    # Create assembly
    my_model.rootAssembly.DatumCsysByDefault(CARTESIAN)
    my_instance = my_model.rootAssembly.Instance(dependent=ON, name=val.INST_NAME, part=my_part)
    my_model.rootAssembly.regenerate()

    my_part.Set(faces=my_part.faces.findAt(((0.5, 0.5, 0.0),)), 
                name=val.BOT_SET_NAME)
    my_part.Set(faces=my_part.faces.findAt(((0.5, 0.5, 1.0),)), 
                name=val.TOP_SET_NAME)
                
    my_part.Set(edges=my_part.edges.findAt(((0.5, 0.0, 0.0),)),
                name=val.XAX_SET_NAME)
    my_part.Set(edges=my_part.edges.findAt(((0.0, 0.5, 0.0),)),
                name=val.YAX_SET_NAME)
    my_part.Set(edges=my_part.edges.findAt(((0.0, 0.0, 0.5),)),
                name=val.ZAX_SET_NAME)
    
    
    # Get corners
    get_corners(my_part)
                
    # Create loading steps
    create_loading(my_model, my_instance)

    # Setup history and field output requests
    my_model.HistoryOutputRequest(
        name='AXIAL_OUTPUT', createStepName='STEP-1',
        region=my_instance.sets['ALL'],
        variables=['U1', 'U2', 'U3', 'RF1', 'RF2', 'RF3'])
        
    my_model.FieldOutputRequest(name='DefaultFieldOutput', 
                                 createStepName='STEP-1')

    # Create job
    my_job = mdb.Job(name=val.JOB_NAME, model=my_model, userSubroutine=val.umatpath)
    
    return my_job


def create_material(the_model):
    the_material = the_model.Material(name=val.MTRL_NAME)
    if val.mtrl_model == 'Chaboche_builtin':
        E       = val.mpar[0]
        nu      = val.mpar[1]
        sy0     = val.mpar[2]
        Hiso    = val.mpar[3]
        kinfinv = val.mpar[4]
        Hkin    = val.mpar[5]
        binfinv = val.mpar[6]
        
        the_material.Elastic(table=((E, nu),))
        the_material.Plastic(table=((sy0, Hkin, Hkin*binfinv),), 
                             hardening=COMBINED, dataType=PARAMETERS, 
                             numBackstresses=1)
        the_material.plastic.CyclicHardening(
                             table=((sy0, 1.0/kinfinv, Hiso*kinfinv),),
                             parameters=ON)
    else:
        the_material.UserMaterial(type=MECHANICAL, unsymm=OFF, 
                                  mechanicalConstants=val.mpar)
        the_material.Depvar(n=val.nstatv)


def create_mesh(the_part, the_cell):
    et = ElemType(elemCode=C3D8, elemLibrary=STANDARD,
            secondOrderAccuracy=OFF, distortionControl=DEFAULT)
    the_part.setElementType(elemTypes=(et,), regions=(the_cell,))
    the_part.seedPart(size=1.0, deviationFactor=0.1, minSizeFactor=0.1)
    the_part.generateMesh()


def get_corners(the_part):
    names = ['ALL', 'F11_PT', 'F22_PT', 'F33_PT', 'F13_PT', 'F23_PT', 'F11+F13_PT', 'F22+F23_PT', 'BOTZ',
                    'P11_PT', 'P22_PT', 'P33_PT']
    o = (0,0,0)
    a = (0,0,1)
    b = (1,0,1)
    c = (0,1,1)
    d = (1,1,1)
    e = (1,0,0)
    f = (0,1,0)
    g = (1,1,0)
    pts = ((o,a,b,c,d,e,f,g),(e,g), (f,g), (a,b,c,d), (a,c), (a,b), (b,d), (c,d), (e,f,g),
            (b,d,e,g), (c,d,f,g), (a,b,c,d))
    the_part.Set(vertices=the_part.vertices.findAt((o,)),name='ORIGIN')
    the_part.Set(vertices=the_part.vertices.findAt((f,)),name='fx0')
    the_part.Set(vertices=the_part.vertices.findAt((e,)),name='ey0')
    for lp,n in zip(pts,names):
        v = []
        for p in lp:
            v.append(the_part.vertices.findAt((p,)))
        
        the_part.Set(vertices=v,name=n)


def create_loading(the_model, the_instance):
    # Fixed points
    # Origin should always be fixed
    the_model.DisplacementBC(name='ORIGIN', createStepName='Initial', 
                             region=the_instance.sets['ORIGIN'], 
                             u1=0.0, u2=0.0, u3=0.0)
                             
    the_model.DisplacementBC(name='BOTZ', createStepName='Initial',
                             region=the_instance.sets['BOTZ'], u3=0.0)
                             
    the_model.DisplacementBC(name='fx0', createStepName='Initial',
                             region=the_instance.sets['fx0'], u1=0.0)
    
    the_model.DisplacementBC(name='ey0', createStepName='Initial',
                             region=the_instance.sets['ey0'], u2=0.0)    
    
    # Load step
    the_model.StaticStep(
        name='STEP-1', previous='Initial',timePeriod=val.total_time,
        timeIncrementationMethod=FIXED, initialInc=val.dt,
        maxNumInc=int(val.total_time/val.dt + 2), nlgeom=val.nlgeom)
    
    # Prescribed points
    the_model.DisplacementBC(name='F11', createStepName='STEP-1',
                             region=the_instance.sets['F11_PT'], 
                             u1=val.F11-1.0)
                             
    the_model.DisplacementBC(name='F22', createStepName='STEP-1',
                             region=the_instance.sets['F22_PT'], 
                             u2=val.F22-1.0)
    
    the_model.DisplacementBC(name='F33', createStepName='STEP-1',
                             region=the_instance.sets['F33_PT'], 
                             u3=val.F33-1.0)
                             
    the_model.DisplacementBC(name='F13', createStepName='STEP-1',
                             region=the_instance.sets['F13_PT'], 
                             u1=val.F13)
                             
    the_model.DisplacementBC(name='F23', createStepName='STEP-1',
                             region=the_instance.sets['F23_PT'], 
                             u2=val.F23)
                             
    the_model.DisplacementBC(name='F11+F13', createStepName='STEP-1',
                             region=the_instance.sets['F11+F13_PT'], 
                             u1=val.F11-1.0+val.F13)
                             
    the_model.DisplacementBC(name='F22+F23', createStepName='STEP-1',
                             region=the_instance.sets['F22+F23_PT'], 
                             u2=val.F22-1.0+val.F23)


                                 
if __name__ == '__main__':              
    setup_model()                              # run the main function