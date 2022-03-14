//
// Created by alban on 11/03/2022.
//

#include "DummyController.h"
#include <sofa/core/ObjectFactory.h>

int DummyControllerClass = sofa::core::RegisterObject("DummyController").add< DummyController >();