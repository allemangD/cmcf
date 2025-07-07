/*==============================================================================

  Program: 3D Slicer

  Portions (c) Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

==============================================================================*/

// .NAME vtkSlicerCMCFlibLogic - slicer logic class for volumes manipulation
// .SECTION Description
// This class manages the logic associated with reading, saving,
// and changing propertied of the volumes

#ifndef __vtkSlicerCMCFlibLogic_h
#define __vtkSlicerCMCFlibLogic_h

// Slicer includes
#include "vtkSlicerModuleLogic.h"

// MRML includes
#include <vtkMRMLModelNode.h>
#include <vtkMRMLScene.h>
#include <vtkMRMLSequenceNode.h>

// STD includes
#include <cstdlib>

#include "vtkSlicerCMCFlibModuleLogicExport.h"

class VTK_SLICER_CMCFLIB_MODULE_LOGIC_EXPORT vtkSlicerCMCFlibLogic : public vtkSlicerModuleLogic {
public:
  static vtkSlicerCMCFlibLogic *New();
  vtkTypeMacro(vtkSlicerCMCFlibLogic, vtkSlicerModuleLogic);
  void PrintSelf(ostream &os, vtkIndent indent) override;

  static void GenerateCMCFSequence(
    vtkMRMLModelNode *model,
    vtkMRMLSequenceNode *sequence,
    double rate = 0.05,
    int stages = 100
  );

  static void IdentifyParabolics(vtkMRMLSequenceNode *sequence, int skip = 5, double tolerance = 0.005);

protected:
  vtkSlicerCMCFlibLogic();
  ~vtkSlicerCMCFlibLogic() override;

  void SetMRMLSceneInternal(vtkMRMLScene *newScene) override;
  /// Register MRML Node classes to Scene. Gets called automatically when the
  /// MRMLScene is attached to this logic class.
  void RegisterNodes() override;
  void UpdateFromMRMLScene() override;
  void OnMRMLSceneNodeAdded(vtkMRMLNode *node) override;
  void OnMRMLSceneNodeRemoved(vtkMRMLNode *node) override;

private:
  vtkSlicerCMCFlibLogic(vtkSlicerCMCFlibLogic const &);  // Not implemented
  void operator=(vtkSlicerCMCFlibLogic const &);         // Not implemented
};

#endif
