#!/usr/bin/env python

#> \file
#> \author Chris Bradley
#> \brief This is an example script to solve a finite elasticity cantilever problem with a growth and constituative law in CellML. The growth occurs just at the bottom of the cantilever in order to cause upward bending. 
#>
#> \section LICENSE
#>
#> Version: MPL 1.1/GPL 2.0/LGPL 2.1
#>
#> The contents of this file are subject to the Mozilla Public License
#> Version 1.1 (the "License"); you may not use this file except in
#> compliance with the License. You may obtain a copy of the License at
#> http://www.mozilla.org/MPL/
#>
#> Software distributed under the License is distributed on an "AS IS"
#> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
#> License for the specific language governing rights and limitations
#> under the License.
#>
#> The Original Code is OpenCMISS
#>
#> The Initial Developer of the Original Code is University of Auckland,
#> Auckland, New Zealand and University of Oxford, Oxford, United
#> Kingdom. Portions created by the University of Auckland and University
#> of Oxford are Copyright (C) 2007 by the University of Auckland and
#> the University of Oxford. All Rights Reserved.
#>
#> Contributor(s): 
#>
#> Alternatively, the contents of this file may be used under the terms of
#> either the GNU General Public License Version 2 or later (the "GPL"), or
#> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
#> in which case the provisions of the GPL or the LGPL are applicable instead
#> of those above. if you wish to allow use of your version of this file only
#> under the terms of either the GPL or the LGPL, and not to allow others to
#> use your version of this file under the terms of the MPL, indicate your
#> decision by deleting the provisions above and replace them with the notice
#> and other provisions required by the GPL or the LGPL. if you do not delete
#> the provisions above, a recipient may use your version of this file under
#> the terms of any one of the MPL, the GPL or the LGPL.
#>

#> Main script
# Add Python bindings directory to PATH
import sys, os

# Intialise OpenCMISS
from opencmiss.opencmiss import OpenCMISS_Python as oc

CONSTANT_LAGRANGE = 0
LINEAR_LAGRANGE = 1
QUADRATIC_LAGRANGE = 2
CUBIC_LAGRANGE = 3

# Set the physical size of the cantilever
width = 10.0
length = 30.0
height = 10.0

# Set the number of elements in the cantilever
numberOfGlobalXElements = 1
numberOfGlobalYElements = 1
numberOfGlobalZElements = 3

# Set the interpolation
uInterpolation = QUADRATIC_LAGRANGE
pInterpolation = LINEAR_LAGRANGE

# Set the growth rates
fibreRate = 0.001 #Or x direction growth
sheetRate = 0.1 #Or y direction growth
normalRate = 0.05 #Or z direction growth

# Set the similation times.
startTime = 0.0
stopTime = 3.0
timeIncrement = 1.0

# materials parameters
c1 = 2.0
c2 = 6.0

force = -0.3

pInit = -6.0
pRef = 0.0

# Set the user numbers
contextUserNumber = 1
coordinateSystemUserNumber = 1
regionUserNumber = 1
uBasisUserNumber = 1
pBasisUserNumber = 2
generatedMeshUserNumber = 1
meshUserNumber = 1
decompositionUserNumber = 1
decomposerUserNumber = 1
geometricFieldUserNumber = 1
fibreFieldUserNumber = 2
dependentFieldUserNumber = 3
equationsSetUserNumber = 1
equationsSetFieldUserNumber = 5
growthCellMLUserNumber = 1
growthCellMLModelsFieldUserNumber = 6
growthCellMLStateFieldUserNumber = 7
growthCellMLParametersFieldUserNumber = 8
constituativeCellMLUserNumber = 2
constituativeCellMLModelsFieldUserNumber = 9
constituativeCellMLParametersFieldUserNumber = 10
constituativeCellMLIntermediateFieldUserNumber = 11
problemUserNumber = 1

numberOfDimensions = 3

if (uInterpolation == LINEAR_LAGRANGE):
    numberOfNodesXi = 2
    numberOfGaussXi = 2
elif (uInterpolation == QUADRATIC_LAGRANGE):
    numberOfNodesXi = 3
    numberOfGaussXi = 3
elif (uInterpolation == CUBIC_LAGRANGE):
    numberOfNodesXi = 4
    numberOfGaussXi = 3
else:
    print('Invalid u interpolation')
    exit()

numberOfXNodes = numberOfGlobalXElements*(numberOfNodesXi-1)+1
numberOfYNodes = numberOfGlobalYElements*(numberOfNodesXi-1)+1
numberOfZNodes = numberOfGlobalZElements*(numberOfNodesXi-1)+1
numberOfNodes = numberOfXNodes*numberOfYNodes*numberOfZNodes

context = oc.Context()
context.Create(contextUserNumber)

worldRegion = oc.Region()
context.WorldRegionGet(worldRegion)

#oc.DiagnosticsSetOn(oc.DiagnosticTypes.FROM,[1,2,3,4,5],"diagnostics",["FiniteElasticity_FiniteElementResidualEvaluate"])

# Get the number of computational nodes and this computational node number
computationEnvironment = oc.ComputationEnvironment()
context.ComputationEnvironmentGet(computationEnvironment)

worldWorkGroup = oc.WorkGroup()
computationEnvironment.WorldWorkGroupGet(worldWorkGroup)
numberOfComputationalNodes = worldWorkGroup.NumberOfGroupNodesGet()
computationalNodeNumber = worldWorkGroup.GroupNodeNumberGet()

# Create a 3D rectangular cartesian coordinate system
coordinateSystem = oc.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber,context)
coordinateSystem.DimensionSet(numberOfDimensions)
coordinateSystem.CreateFinish()

# Create a region and assign the coordinate system to the region
region = oc.Region()
region.CreateStart(regionUserNumber,worldRegion)
region.LabelSet("Region")
region.CoordinateSystemSet(coordinateSystem)
region.CreateFinish()

# Define basis functions

uBasis = oc.Basis()
uBasis.CreateStart(uBasisUserNumber,context)
uBasis.NumberOfXiSet(numberOfDimensions)
uBasis.TypeSet(oc.BasisTypes.LAGRANGE_HERMITE_TP)
if (uInterpolation == LINEAR_LAGRANGE):
    uBasis.InterpolationXiSet([oc.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*numberOfDimensions)
elif (uInterpolation == QUADRATIC_LAGRANGE):
    uBasis.InterpolationXiSet([oc.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE]*numberOfDimensions)
elif (uInterpolation == CUBIC_LAGRANGE):
    uBasis.InterpolationXiSet([oc.BasisInterpolationSpecifications.CUBIC_LAGRANGE]*numberOfDimensions)
else:
    print('Invalid u interpolation')
    exit()
uBasis.QuadratureNumberOfGaussXiSet([numberOfGaussXi]*numberOfDimensions)
uBasis.CreateFinish()

if (pInterpolation > CONSTANT_LAGRANGE):
    pBasis = oc.Basis()
    pBasis.CreateStart(pBasisUserNumber,context)
    pBasis.NumberOfXiSet(numberOfDimensions)
    pBasis.TypeSet(oc.BasisTypes.LAGRANGE_HERMITE_TP)
    if (pInterpolation == LINEAR_LAGRANGE):
        pBasis.InterpolationXiSet([oc.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*numberOfDimensions)
    elif (pInterpolation == QUADRATIC_LAGRANGE):
        pBasis.InterpolationXiSet([oc.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE]*numberOfDimensions)
    else:
        print('Invalid p interpolation')
        exit()
    pBasis.QuadratureNumberOfGaussXiSet([numberOfGaussXi]*numberOfDimensions)
    pBasis.CreateFinish()

# Start the creation of a generated mesh in the region
generatedMesh = oc.GeneratedMesh()
generatedMesh.CreateStart(generatedMeshUserNumber,region)
generatedMesh.TypeSet(oc.GeneratedMeshTypes.REGULAR)
if (pInterpolation == CONSTANT_LAGRANGE):
    generatedMesh.BasisSet([uBasis])
else:
    generatedMesh.BasisSet([uBasis,pBasis])
generatedMesh.ExtentSet([width,height,length])
generatedMesh.NumberOfElementsSet([numberOfGlobalXElements,numberOfGlobalYElements,numberOfGlobalZElements])
# Finish the creation of a generated mesh in the region
mesh = oc.Mesh()
generatedMesh.CreateFinish(meshUserNumber,mesh)

# Create a decomposition for the mesh
decomposition = oc.Decomposition()
decomposition.CreateStart(decompositionUserNumber,mesh)
decomposition.CreateFinish()

# Decompose 
decomposer = oc.Decomposer()
decomposer.CreateStart(decomposerUserNumber,worldRegion,worldWorkGroup)
decompositionIndex = decomposer.DecompositionAdd(decomposition)
decomposer.CreateFinish()

# Create a field for the geometry
geometricField = oc.Field()
geometricField.CreateStart(geometricFieldUserNumber,region)
geometricField.DecompositionSet(decomposition)
geometricField.TypeSet(oc.FieldTypes.GEOMETRIC)
geometricField.VariableLabelSet(oc.FieldVariableTypes.U,"Geometry")
geometricField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,1,1)
geometricField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,2,1)
geometricField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,3,1)
geometricField.ScalingTypeSet(oc.FieldScalingTypes.ARITHMETIC_MEAN)
geometricField.CreateFinish()

# Update the geometric field parameters from generated mesh
generatedMesh.GeometricParametersCalculate(geometricField)

# Create a fibre field and attach it to the geometric field
fibreField = oc.Field()
fibreField.CreateStart(fibreFieldUserNumber,region)
fibreField.TypeSet(oc.FieldTypes.FIBRE)
fibreField.DecompositionSet(decomposition)
fibreField.GeometricFieldSet(geometricField)
fibreField.VariableLabelSet(oc.FieldVariableTypes.U,"Fibre")
fibreField.ScalingTypeSet(oc.FieldScalingTypes.ARITHMETIC_MEAN)
fibreField.CreateFinish()

# Create the dependent field
dependentField = oc.Field()
dependentField.CreateStart(dependentFieldUserNumber,region)
dependentField.TypeSet(oc.FieldTypes.GEOMETRIC_GENERAL)  
dependentField.DecompositionSet(decomposition)
dependentField.GeometricFieldSet(geometricField) 
dependentField.DependentTypeSet(oc.FieldDependentTypes.DEPENDENT) 
# Set the field to have 5 variables: U - dependent; del U/del n - tractions; U1 - strain; U2 - stress; U3 - growth
dependentField.NumberOfVariablesSet(5)
dependentField.VariableTypesSet([oc.FieldVariableTypes.U,oc.FieldVariableTypes.T,oc.FieldVariableTypes.U1,oc.FieldVariableTypes.U2,oc.FieldVariableTypes.U3])
dependentField.VariableLabelSet(oc.FieldVariableTypes.U,"Displacement")
dependentField.VariableLabelSet(oc.FieldVariableTypes.T,"Traction")
dependentField.VariableLabelSet(oc.FieldVariableTypes.U1,"Strain")
dependentField.VariableLabelSet(oc.FieldVariableTypes.U2,"Stress")
dependentField.VariableLabelSet(oc.FieldVariableTypes.U3,"Growth")
dependentField.NumberOfComponentsSet(oc.FieldVariableTypes.U,4)
dependentField.NumberOfComponentsSet(oc.FieldVariableTypes.T,4)
dependentField.NumberOfComponentsSet(oc.FieldVariableTypes.U1,6)
dependentField.NumberOfComponentsSet(oc.FieldVariableTypes.U2,6)
dependentField.NumberOfComponentsSet(oc.FieldVariableTypes.U3,3)
if (pInterpolation == CONSTANT_LAGRANGE):
    dependentField.ComponentInterpolationSet(oc.FieldVariableTypes.U,4,oc.FieldInterpolationTypes.ELEMENT_BASED)
    dependentField.ComponentInterpolationSet(oc.FieldVariableTypes.T,4,oc.FieldInterpolationTypes.ELEMENT_BASED)
dependentField.ComponentInterpolationSet(oc.FieldVariableTypes.U1,1,oc.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(oc.FieldVariableTypes.U1,2,oc.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(oc.FieldVariableTypes.U1,3,oc.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(oc.FieldVariableTypes.U1,4,oc.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(oc.FieldVariableTypes.U1,5,oc.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(oc.FieldVariableTypes.U1,6,oc.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(oc.FieldVariableTypes.U2,1,oc.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(oc.FieldVariableTypes.U2,2,oc.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(oc.FieldVariableTypes.U2,3,oc.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(oc.FieldVariableTypes.U2,4,oc.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(oc.FieldVariableTypes.U2,5,oc.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(oc.FieldVariableTypes.U2,6,oc.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(oc.FieldVariableTypes.U3,1,oc.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(oc.FieldVariableTypes.U3,2,oc.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(oc.FieldVariableTypes.U3,3,oc.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ScalingTypeSet(oc.FieldScalingTypes.ARITHMETIC_MEAN)
dependentField.CreateFinish()

# Initialise dependent field from undeformed geometry
oc.Field.ParametersToFieldParametersComponentCopy(
    geometricField,oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,1,
    dependentField,oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,1)
oc.Field.ParametersToFieldParametersComponentCopy(
    geometricField,oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,2,
    dependentField,oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,2)
oc.Field.ParametersToFieldParametersComponentCopy(
    geometricField,oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,3,
    dependentField,oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,3)
# Initialise the hydrostatic pressure
oc.Field.ComponentValuesInitialiseDP(
    dependentField,oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,4,pInit)

# Update the dependent field
dependentField.ParameterSetUpdateStart(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES)
dependentField.ParameterSetUpdateFinish(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES)

# Create the equations_set
equationsSetField = oc.Field()
equationsSet = oc.EquationsSet()
equationsSetSpecification = [oc.EquationsSetClasses.ELASTICITY,
    oc.EquationsSetTypes.FINITE_ELASTICITY,
    oc.EquationsSetSubtypes.CONSTIT_AND_GROWTH_LAW_IN_CELLML]
equationsSet.CreateStart(equationsSetUserNumber,region,fibreField,
                         equationsSetSpecification,equationsSetFieldUserNumber,equationsSetField)
equationsSet.CreateFinish()

equationsSet.DependentCreateStart(dependentFieldUserNumber,dependentField)
equationsSet.DependentCreateFinish()

# Create the CellML environment for the growth law. Set the rates as known so that we can spatially vary them.
growthCellML = oc.CellML()
growthCellML.CreateStart(growthCellMLUserNumber,region)
growthCellMLIdx = growthCellML.ModelImport("stressgrowth.cellml")
growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/bff")
growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/bss")
growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/bnn")
growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/S11")
growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/S22")
growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/S33")
growthCellML.CreateFinish()

# Create CellML <--> OpenCMISS field maps. Map the lambda's to the U3/growth dependent field variable
growthCellML.FieldMapsCreateStart()
growthCellML.CreateFieldToCellMLMap(dependentField,oc.FieldVariableTypes.U2,1,oc.FieldParameterSetTypes.VALUES,
	                            growthCellMLIdx,"Main/S11",oc.FieldParameterSetTypes.VALUES)
growthCellML.CreateFieldToCellMLMap(dependentField,oc.FieldVariableTypes.U2,2,oc.FieldParameterSetTypes.VALUES,
	                            growthCellMLIdx,"Main/S22",oc.FieldParameterSetTypes.VALUES)
growthCellML.CreateFieldToCellMLMap(dependentField,oc.FieldVariableTypes.U2,3,oc.FieldParameterSetTypes.VALUES,
	                            growthCellMLIdx,"Main/S33",oc.FieldParameterSetTypes.VALUES)
growthCellML.CreateCellMLToFieldMap(growthCellMLIdx,"Main/lambda1",oc.FieldParameterSetTypes.VALUES,
                                    dependentField,oc.FieldVariableTypes.U3,1,oc.FieldParameterSetTypes.VALUES)
growthCellML.CreateCellMLToFieldMap(growthCellMLIdx,"Main/lambda2",oc.FieldParameterSetTypes.VALUES,
                                    dependentField,oc.FieldVariableTypes.U3,2,oc.FieldParameterSetTypes.VALUES)
growthCellML.CreateCellMLToFieldMap(growthCellMLIdx,"Main/lambda3",oc.FieldParameterSetTypes.VALUES,
                                    dependentField,oc.FieldVariableTypes.U3,3,oc.FieldParameterSetTypes.VALUES)
growthCellML.FieldMapsCreateFinish()

# Create the CELL models field
growthCellMLModelsField = oc.Field()
growthCellML.ModelsFieldCreateStart(growthCellMLModelsFieldUserNumber,growthCellMLModelsField)
growthCellMLModelsField.VariableLabelSet(oc.FieldVariableTypes.U,"GrowthModelMap")
growthCellML.ModelsFieldCreateFinish()

# Create the CELL parameters field
growthCellMLParametersField = oc.Field()
growthCellML.ParametersFieldCreateStart(growthCellMLParametersFieldUserNumber,growthCellMLParametersField)
growthCellMLParametersField.VariableLabelSet(oc.FieldVariableTypes.U,"GrowthParameters")
growthCellML.ParametersFieldCreateFinish()

# Set the growth rates
fibreRateComponentNumber = growthCellML.FieldComponentGet(growthCellMLIdx,oc.CellMLFieldTypes.PARAMETERS,"Main/bff")
sheetRateComponentNumber = growthCellML.FieldComponentGet(growthCellMLIdx,oc.CellMLFieldTypes.PARAMETERS,"Main/bss")
normalRateComponentNumber = growthCellML.FieldComponentGet(growthCellMLIdx,oc.CellMLFieldTypes.PARAMETERS,"Main/bnn")
growthCellMLParametersField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                                        fibreRateComponentNumber,fibreRate)
growthCellMLParametersField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                                        sheetRateComponentNumber,sheetRate)
growthCellMLParametersField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                                        normalRateComponentNumber,normalRate)

# Create the CELL state field
growthCellMLStateField = oc.Field()
growthCellML.StateFieldCreateStart(growthCellMLStateFieldUserNumber,growthCellMLStateField)
growthCellMLStateField.VariableLabelSet(oc.FieldVariableTypes.U,"GrowthState")
growthCellML.StateFieldCreateFinish()

# Create the CellML environment for the consitutative law
constituativeCellML = oc.CellML()
constituativeCellML.CreateStart(constituativeCellMLUserNumber,region)
constituativeCellMLIdx = constituativeCellML.ModelImport("mooneyrivlin.cellml")
constituativeCellML.VariableSetAsKnown(constituativeCellMLIdx,"equations/C11")
constituativeCellML.VariableSetAsKnown(constituativeCellMLIdx,"equations/C12")
constituativeCellML.VariableSetAsKnown(constituativeCellMLIdx,"equations/C13")
constituativeCellML.VariableSetAsKnown(constituativeCellMLIdx,"equations/C22")
constituativeCellML.VariableSetAsKnown(constituativeCellMLIdx,"equations/C23")
constituativeCellML.VariableSetAsKnown(constituativeCellMLIdx,"equations/C33")
constituativeCellML.VariableSetAsKnown(constituativeCellMLIdx,"equations/c1")
constituativeCellML.VariableSetAsKnown(constituativeCellMLIdx,"equations/c2")
constituativeCellML.VariableSetAsWanted(constituativeCellMLIdx,"equations/Tdev11")
constituativeCellML.VariableSetAsWanted(constituativeCellMLIdx,"equations/Tdev12")
constituativeCellML.VariableSetAsWanted(constituativeCellMLIdx,"equations/Tdev13")
constituativeCellML.VariableSetAsWanted(constituativeCellMLIdx,"equations/Tdev22")
constituativeCellML.VariableSetAsWanted(constituativeCellMLIdx,"equations/Tdev23")
constituativeCellML.VariableSetAsWanted(constituativeCellMLIdx,"equations/Tdev33")
constituativeCellML.CreateFinish()

# Create CellML <--> OpenCMISS field maps. Map the stress and strain fields.
constituativeCellML.FieldMapsCreateStart()
constituativeCellML.CreateFieldToCellMLMap(dependentField,oc.FieldVariableTypes.U1,1,oc.FieldParameterSetTypes.VALUES,
    constituativeCellMLIdx,"equations/C11",oc.FieldParameterSetTypes.VALUES)
constituativeCellML.CreateFieldToCellMLMap(dependentField,oc.FieldVariableTypes.U1,2,oc.FieldParameterSetTypes.VALUES,
    constituativeCellMLIdx,"equations/C12",oc.FieldParameterSetTypes.VALUES)
constituativeCellML.CreateFieldToCellMLMap(dependentField,oc.FieldVariableTypes.U1,3,oc.FieldParameterSetTypes.VALUES,
    constituativeCellMLIdx,"equations/C13",oc.FieldParameterSetTypes.VALUES)
constituativeCellML.CreateFieldToCellMLMap(dependentField,oc.FieldVariableTypes.U1,4,oc.FieldParameterSetTypes.VALUES,
    constituativeCellMLIdx,"equations/C22",oc.FieldParameterSetTypes.VALUES)
constituativeCellML.CreateFieldToCellMLMap(dependentField,oc.FieldVariableTypes.U1,5,oc.FieldParameterSetTypes.VALUES,
    constituativeCellMLIdx,"equations/C23",oc.FieldParameterSetTypes.VALUES)
constituativeCellML.CreateFieldToCellMLMap(dependentField,oc.FieldVariableTypes.U1,6,oc.FieldParameterSetTypes.VALUES,
    constituativeCellMLIdx,"equations/C33",oc.FieldParameterSetTypes.VALUES)
constituativeCellML.CreateCellMLToFieldMap(constituativeCellMLIdx,"equations/Tdev11",oc.FieldParameterSetTypes.VALUES,
    dependentField,oc.FieldVariableTypes.U2,1,oc.FieldParameterSetTypes.VALUES)
constituativeCellML.CreateCellMLToFieldMap(constituativeCellMLIdx,"equations/Tdev12",oc.FieldParameterSetTypes.VALUES,
    dependentField,oc.FieldVariableTypes.U2,2,oc.FieldParameterSetTypes.VALUES)
constituativeCellML.CreateCellMLToFieldMap(constituativeCellMLIdx,"equations/Tdev13",oc.FieldParameterSetTypes.VALUES,
    dependentField,oc.FieldVariableTypes.U2,3,oc.FieldParameterSetTypes.VALUES)
constituativeCellML.CreateCellMLToFieldMap(constituativeCellMLIdx,"equations/Tdev22",oc.FieldParameterSetTypes.VALUES,
    dependentField,oc.FieldVariableTypes.U2,4,oc.FieldParameterSetTypes.VALUES)
constituativeCellML.CreateCellMLToFieldMap(constituativeCellMLIdx,"equations/Tdev23",oc.FieldParameterSetTypes.VALUES,
    dependentField,oc.FieldVariableTypes.U2,5,oc.FieldParameterSetTypes.VALUES)
constituativeCellML.CreateCellMLToFieldMap(constituativeCellMLIdx,"equations/Tdev33",oc.FieldParameterSetTypes.VALUES,
    dependentField,oc.FieldVariableTypes.U2,6,oc.FieldParameterSetTypes.VALUES)
constituativeCellML.FieldMapsCreateFinish()

# Create the CELL models field
constituativeCellMLModelsField = oc.Field()
constituativeCellML.ModelsFieldCreateStart(constituativeCellMLModelsFieldUserNumber,constituativeCellMLModelsField)
constituativeCellMLModelsField.VariableLabelSet(oc.FieldVariableTypes.U,"ConstituativeModelMap")
constituativeCellML.ModelsFieldCreateFinish()

# Create the CELL parameters field
constituativeCellMLParametersField = oc.Field()
constituativeCellML.ParametersFieldCreateStart(constituativeCellMLParametersFieldUserNumber,constituativeCellMLParametersField)
constituativeCellMLParametersField.VariableLabelSet(oc.FieldVariableTypes.U,"ConstituativeParameters")
constituativeCellML.ParametersFieldCreateFinish()

# Set up the materials constants
c1ComponentNumber = constituativeCellML.FieldComponentGet(constituativeCellMLIdx,oc.CellMLFieldTypes.PARAMETERS,"equations/c1")
c2ComponentNumber = constituativeCellML.FieldComponentGet(constituativeCellMLIdx,oc.CellMLFieldTypes.PARAMETERS,"equations/c2")
constituativeCellMLParametersField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                                               c1ComponentNumber,c1)
constituativeCellMLParametersField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                                               c2ComponentNumber,c2)

# Create the CELL intermediate field
constituativeCellMLIntermediateField = oc.Field()
constituativeCellML.IntermediateFieldCreateStart(constituativeCellMLIntermediateFieldUserNumber,constituativeCellMLIntermediateField)
constituativeCellMLIntermediateField.VariableLabelSet(oc.FieldVariableTypes.U,"ConstituativeIntermediate")
constituativeCellML.IntermediateFieldCreateFinish()

# Create equations
equations = oc.Equations()
equationsSet.EquationsCreateStart(equations)
equations.sparsityType = oc.EquationsSparsityTypes.SPARSE
equations.outputType = oc.EquationsOutputTypes.NONE
equationsSet.EquationsCreateFinish()

# Define the problem
problem = oc.Problem()
problemSpecification = [oc.ProblemClasses.ELASTICITY,
        oc.ProblemTypes.FINITE_ELASTICITY,
        oc.ProblemSubtypes.FINITE_ELASTICITY_WITH_GROWTH_CELLML]
problem.CreateStart(problemUserNumber,context,problemSpecification)
problem.CreateFinish()

# Create control loops
timeLoop = oc.ControlLoop()
problem.ControlLoopCreateStart()
problem.ControlLoopGet([oc.ControlLoopIdentifiers.NODE],timeLoop)
timeLoop.TimesSet(startTime,stopTime,timeIncrement)
problem.ControlLoopCreateFinish()

# Create problem solvers
odeIntegrationSolver = oc.Solver()
nonlinearSolver = oc.Solver()
linearSolver = oc.Solver()
cellMLEvaluationSolver = oc.Solver()
problem.SolversCreateStart()
problem.SolverGet([oc.ControlLoopIdentifiers.NODE],1,odeIntegrationSolver)
problem.SolverGet([oc.ControlLoopIdentifiers.NODE],2,nonlinearSolver)
nonlinearSolver.outputType = oc.SolverOutputTypes.MONITOR
nonlinearSolver.NewtonJacobianCalculationTypeSet(oc.JacobianCalculationTypes.FD)
nonlinearSolver.NewtonAbsoluteToleranceSet(1e-11)
nonlinearSolver.NewtonSolutionToleranceSet(1e-11)
nonlinearSolver.NewtonRelativeToleranceSet(1e-11)
nonlinearSolver.NewtonCellMLSolverGet(cellMLEvaluationSolver)
nonlinearSolver.NewtonLinearSolverGet(linearSolver)
linearSolver.linearType = oc.LinearSolverTypes.DIRECT
problem.SolversCreateFinish()

# Create nonlinear equations and add equations set to solver equations
nonlinearEquations = oc.SolverEquations()
problem.SolverEquationsCreateStart()
nonlinearSolver.SolverEquationsGet(nonlinearEquations)
nonlinearEquations.sparsityType = oc.SolverEquationsSparsityTypes.SPARSE
nonlinearEquationsSetIndex = nonlinearEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

# Create CellML equations and add growth and constituative equations to the solvers
growthEquations = oc.CellMLEquations()
constituativeEquations = oc.CellMLEquations()
problem.CellMLEquationsCreateStart()
odeIntegrationSolver.CellMLEquationsGet(growthEquations)
growthEquationsIndex = growthEquations.CellMLAdd(growthCellML)
cellMLEvaluationSolver.CellMLEquationsGet(constituativeEquations)
constituativeEquationsIndex = constituativeEquations.CellMLAdd(constituativeCellML)
problem.CellMLEquationsCreateFinish()

# Prescribe boundary conditions (absolute nodal parameters)
boundaryConditions = oc.BoundaryConditions()
nonlinearEquations.BoundaryConditionsCreateStart(boundaryConditions)

for widthNodeIdx in range(1,numberOfXNodes+1):
    for heightNodeIdx in range(1,numberOfYNodes+1):
        # Set left hand build in nodes ot no displacement
        nodeIdx=widthNodeIdx+(heightNodeIdx-1)*numberOfXNodes
        boundaryConditions.AddNode(dependentField,oc.FieldVariableTypes.U,1,1,nodeIdx,1,
                                   oc.BoundaryConditionsTypes.FIXED,0.0)
        boundaryConditions.AddNode(dependentField,oc.FieldVariableTypes.U,1,1,nodeIdx,2,
                                   oc.BoundaryConditionsTypes.FIXED,0.0)
        boundaryConditions.AddNode(dependentField,oc.FieldVariableTypes.U,1,1,nodeIdx,3,
                                   oc.BoundaryConditionsTypes.FIXED,0.0)
    # Set downward force on right-hand edge
    nodeIdx=numberOfNodes-widthNodeIdx+1
    boundaryConditions.AddNode(dependentField,oc.FieldVariableTypes.T,1,1,nodeIdx,2,
                               oc.BoundaryConditionsTypes.NEUMANN_POINT,force)
# Set reference pressure
boundaryConditions.AddNode(dependentField,oc.FieldVariableTypes.U,1,1,numberOfNodes,4,
                           oc.BoundaryConditionsTypes.FIXED,pRef)
 
nonlinearEquations.BoundaryConditionsCreateFinish()

# Solve the problem
problem.Solve()

if not os.path.exists("./results"):
    os.makedirs("./results")

# Export results
fields = oc.Fields()
fields.CreateRegion(region)
fields.NodesExport("./results/CantileverGrowth","FORTRAN")
fields.ElementsExport("./results/CantileverGrowth","FORTRAN")
fields.Finalise()

