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

// CMCFlib Logic includes
#include <vtkSlicerCMCFlibLogic.h>

// CMCFlib includes
#include "qSlicerCMCFlibModule.h"

//-----------------------------------------------------------------------------
class qSlicerCMCFlibModulePrivate {
public:
  qSlicerCMCFlibModulePrivate();
};
auto
//-----------------------------------------------------------------------------
// qSlicerCMCFlibModulePrivate methods

//-----------------------------------------------------------------------------
qSlicerCMCFlibModulePrivate::qSlicerCMCFlibModulePrivate() {}

//-----------------------------------------------------------------------------
// qSlicerCMCFlibModule methods

//-----------------------------------------------------------------------------
qSlicerCMCFlibModule::qSlicerCMCFlibModule(QObject *_parent)
  : Superclass(_parent), d_ptr(new qSlicerCMCFlibModulePrivate) {
  this->setWidgetRepresentationCreationEnabled(false);
}

//-----------------------------------------------------------------------------
qSlicerCMCFlibModule::~qSlicerCMCFlibModule() {}

//-----------------------------------------------------------------------------
QString qSlicerCMCFlibModule::helpText() const {
  return "This is a loadable module that can be bundled in an extension";
}

//-----------------------------------------------------------------------------
QString qSlicerCMCFlibModule::acknowledgementText() const {
  return "This work was partially funded by NIH grant NXNNXXNNNNNN-NNXN";
}

//-----------------------------------------------------------------------------
QStringList qSlicerCMCFlibModule::contributors() const {
  QStringList moduleContributors;
  moduleContributors << QString("John Doe (AnyWare Corp.)");
  return moduleContributors;
}

//-----------------------------------------------------------------------------
QIcon qSlicerCMCFlibModule::icon() const { return QIcon(":/Icons/CMCFlib.png"); }

//-----------------------------------------------------------------------------
QStringList qSlicerCMCFlibModule::categories() const { return QStringList() << "Examples"; }

//-----------------------------------------------------------------------------
QStringList qSlicerCMCFlibModule::dependencies() const { return QStringList(); }

//-----------------------------------------------------------------------------
bool qSlicerCMCFlibModule::isHidden() const { return true; }

//-----------------------------------------------------------------------------
void qSlicerCMCFlibModule::setup() { this->Superclass::setup(); }

//-----------------------------------------------------------------------------
qSlicerAbstractModuleRepresentation *
qSlicerCMCFlibModule::createWidgetRepresentation() { return nullptr; }

//-----------------------------------------------------------------------------
vtkMRMLAbstractLogic *qSlicerCMCFlibModule::createLogic() { return vtkSlicerCMCFlibLogic::New(); }