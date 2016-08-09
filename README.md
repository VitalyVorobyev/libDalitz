# libDalitz
Set of tools for Dalitz plot analysis of a threebody particle decay

Introduction

This project is c++ shared library for Dalitz analysis of a three-body particle decay. The following features are implemented at the moment:

    Description of threebody decay phase space in terms of Dalitz variables;
    Isobar model of a decay amplitude containing arbitrary numbers of scalar, vector and tensor resonanses;
    Calculation of parameters values for model-independent analysis with equal-phase binning;
    MC integration of the decay probability over the Dalitz plot;
    Generator of Dalitz plot distribution for a specific decay model.

List of classes

    DalitzPhaseSpace
    DalitzModel
    SymDalitzModel
    KspipiModel
    B0toD0pipiModel
    RandomDalitzPoint
    DalitzMCIntegral
    DalitzGenerator
    ModelIntegral
# libTatami
