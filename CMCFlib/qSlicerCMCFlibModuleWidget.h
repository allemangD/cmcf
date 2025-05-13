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

#ifndef __qSlicerCMCFlibModuleWidget_h
#define __qSlicerCMCFlibModuleWidget_h

// Slicer includes
#include "qSlicerAbstractModuleWidget.h"

#include "qSlicerCMCFlibModuleExport.h"

class qSlicerCMCFlibModuleWidgetPrivate;
class vtkMRMLNode;

class Q_SLICER_QTMODULES_CMCFLIB_EXPORT qSlicerCMCFlibModuleWidget :
  public qSlicerAbstractModuleWidget
{
  Q_OBJECT

public:

  typedef qSlicerAbstractModuleWidget Superclass;
  qSlicerCMCFlibModuleWidget(QWidget *parent=0);
  virtual ~qSlicerCMCFlibModuleWidget();

public slots:


protected:
  QScopedPointer<qSlicerCMCFlibModuleWidgetPrivate> d_ptr;

  void setup() override;

private:
  Q_DECLARE_PRIVATE(qSlicerCMCFlibModuleWidget);
  Q_DISABLE_COPY(qSlicerCMCFlibModuleWidget);
};

#endif
