from sketch         import *
from part           import *
from material       import *
from section        import *
from assembly       import *
from step           import *
from load           import *
from mesh           import *
from job            import *

import mmf_abaver_model_val as val

def setup_model():
    my_model = mdb.models.values()[0]

    # Create sketch 
    my_sketch = create_sketch(my_model)

    # Create the solid/shell
    my_part = my_model.Part(dimensionality=AXISYMMETRIC, name=val.PART_NAME,
                            type=DEFORMABLE_BODY, twist=True)
    my_part.BaseShell(sketch=my_sketch)
    del my_sketch

    # Create the material
    my_material = create_material(my_model)

    # Section
    my_section = my_model.HomogeneousSolidSection(
        name='Section-1', material=val.MTRL_NAME, thickness=None)
    my_face = my_part.faces.findAt((( val.rm, val.hm, 0.),),)
    my_part.SectionAssignment(
        region=Region(faces = my_face), 
        sectionName='Section-1', 
        thicknessAssignment=FROM_GEOMETRY)
        
    # Create sets
    my_part.Set(edges=my_part.edges.findAt(((val.rm, 0.0, 0.0),)), 
                name=val.BOT_SET_NAME)
    my_part.Set(edges=my_part.edges.findAt(((val.rm, val.h0, 0.0),)), 
                name=val.TOP_SET_NAME)
    my_part.Set(edges=my_part.edges.findAt(((val.ri, val.hm, 0.0),)), 
                name=val.INS_SET_NAME)
    my_part.Set(edges=my_part.edges.findAt(((val.ro, val.hm, 0.0),)), 
                name=val.OUT_SET_NAME)
                
    # Create mesh
    create_mesh(my_part, my_face)

    # Create assembly
    my_model.rootAssembly.DatumCsysByDefault(CARTESIAN)
    my_instance = my_model.rootAssembly.Instance(dependent=ON, name=val.INST_NAME, part=my_part)
    my_model.rootAssembly.regenerate()

    # Create loading steps
    create_loading(my_model, my_instance)

    # Setup history and field output requests
    setup_output_requests(my_model, my_instance)

    # Create job
    my_job = mdb.Job(name=val.JOB_NAME, model=my_model, userSubroutine=val.umatpath)
    
    return my_job

def create_sketch(the_model):
    the_sketch = the_model.ConstrainedSketch(
        name='__profile__', sheetSize=20)
        
    #Centerline needed for axisymmetric parts:
    the_sketch.ConstructionLine(point1=(0.0,0.0), 
                                point2=(0.0, val.h0))
                                
    # Define the rectangle shape
    p1 = (val.ri,  0.0)
    p2 = (val.ro, 0.0)
    p3 = (val.ro, val.h0)
    p4 = (val.ri, val.h0)
    
    the_sketch.Line(point1=p1, point2=p2)
    the_sketch.Line(point1=p2, point2=p3)
    the_sketch.Line(point1=p3, point2=p4)
    the_sketch.Line(point1=p4, point2=p1)
    return the_sketch

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
        
def create_mesh(the_part, the_face):
    the_part.setMeshControls(elemShape=QUAD, regions=the_face, technique=STRUCTURED)
    if val.element_order == 1:
        et = ElemType(elemCode=CGAX4, elemLibrary=STANDARD)
    else:
        et = ElemType(elemCode=CGAX8, elemLibrary=STANDARD)

    the_part.setElementType(elemTypes=(et,), regions=(the_face,))

    # val.nel elements in radial direction
    the_part.seedEdgeByNumber(
        edges=the_part.edges.findAt(((val.rm, 0.0, 0.0),)),
        number = val.nel, constraint=FIXED)
    # One element in axial direction
    the_part.seedEdgeByNumber(
        edges=the_part.edges.findAt(((val.ri, val.hm, 0.0),)),
        number = 1, constraint=FIXED)
        
    the_part.generateMesh()

def create_loading(the_model, the_instance):
    # Bottom should be fixed throughout
    the_model.DisplacementBC(name='FIXED_BC', createStepName='Initial', 
                             region=the_instance.sets[val.BOT_SET_NAME], 
                             u2=0.0, ur2=0.0)
    # Create amplitudes
    for i in range(len(val.amp)):
        the_model.TabularAmplitude(
            name='AMP'+str(i+1),data=val.amp[i],smooth=0.0, timeSpan=TOTAL)
    
    # Create steps
    prev_step = 'Initial'
    for i in range(len(val.step_time)):
        the_model.StaticStep(
            name='STEP'+str(i+1), previous=prev_step,
            timePeriod=val.step_time[i], 
            timeIncrementationMethod=FIXED, initialInc=val.dt,
            maxNumInc=int(val.step_time[i]/val.dt + 2),
            nlgeom=val.nlgeom)
        prev_step = 'STEP'+str(i+1)
        
    # Apply prescribed displacements
    the_model.DisplacementBC(name='AXIAL_DISP', createStepName='STEP1', 
                             region=the_instance.sets[val.TOP_SET_NAME],
                             u2=val.epsz_amp*val.h0, 
                             amplitude='AMP1')
    the_model.DisplacementBC(name='ROTATION', createStepName='STEP1',
                             region=the_instance.sets[val.TOP_SET_NAME],
                             ur2=val.gamma_amp*val.h0/val.ro, 
                             amplitude='AMP2')
    the_model.DisplacementBC(name='INSIDE', createStepName='STEP1',
                             region=the_instance.sets[val.INS_SET_NAME],
                             u1=val.cstri_amp*val.ri, 
                             amplitude='AMP3')
    the_model.DisplacementBC(name='OUTSIDE', createStepName='STEP1',
                             region=the_instance.sets[val.OUT_SET_NAME],
                             u1=val.cstro_amp*val.ro, 
                             amplitude='AMP4')
    
def setup_output_requests(the_model, the_instance):
    the_model.HistoryOutputRequest(
        name='AXIAL_OUTPUT', createStepName='STEP1', 
        region=the_instance.sets[val.TOP_SET_NAME], 
        variables=['U2', 'UR2', 'RF2', 'RM2',])
    the_model.HistoryOutputRequest(
        name='INSIDE_OUTPUT', createStepName='STEP1', 
        region=the_instance.sets[val.INS_SET_NAME], 
        variables=['U1', 'RF1',])
    the_model.HistoryOutputRequest(
        name='OUTSIDE_OUTPUT', createStepName='STEP1', 
        region=the_instance.sets[val.OUT_SET_NAME], 
        variables=['U1', 'RF1'])
    the_model.HistoryOutputRequest(
        name='RADIAL_OUTPUT', createStepName='STEP1', 
        region=the_instance.sets[val.TOP_SET_NAME], 
        variables=['U1'])
    the_model.FieldOutputRequest(name='DefaultFieldOutput', 
                                 createStepName='STEP1')
        
if __name__ == '__main__':              
    setup_model()                              # run the main function