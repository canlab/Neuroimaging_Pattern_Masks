% this script converts ants affine transformation matrices to SPM format transformation matrices
% afftransform is the AffineTransform_double_3_3 variable of ants *.mat files
% m_Center is the fixed variable of ants *.mat files
% this script was obtained here: 
% https://github.com/netstim/leaddbs/blob/master/helpers/ea_antsmat2mat.m
% more details are also available here:
% https://github.com/ANTsX/ANTs/wiki/ITK-affine-transform-conversion

function mat=ea_antsmat2mat(afftransform,m_Center)
% followed instructions from
% https://www.neuro.polymtl.ca/tips_and_tricks/how_to_use_ants and ITK code
% in ComputeOffset() itkMatrixOffsetTransformBase.hxx

mat=[reshape(afftransform(1:9),[3,3])',afftransform(10:12)];

m_Translation=mat(:,4);
mat=[mat;[0,0,0,1]];

for i=1:3
    m_Offset(i) = m_Translation(i) + m_Center(i);
    for j=1:3
       m_Offset(i) = m_Offset(i)-(mat(i,j) * m_Center(j));  % (i,j) should remain the same since in C indexing rows before cols, too.
    end
end

mat(1:3,4)=m_Offset;
mat=inv(mat);

% convert RAS to LPS (ITK uses LPS)
mat=mat.*...
    [1  1 -1 -1
     1  1 -1 -1
    -1 -1  1  1
     1  1  1  1];

% original code in itkMatrixOffsetTransformBase > ComputeOffset
% {
%   const MatrixType & matrix = this->GetMatrix();
%
%   OffsetType offset;
%   for(unsigned int i=0; i<NOutputDimensions; i++)
%     {
%     offset[i] = m_Translation[i] + m_Center[i];
%     for(unsigned int j=0; j<NInputDimensions; j++)
%       {
%       offset[i] -= matrix[i][j] * m_Center[j];
%       }
%     }
%
%   m_Offset = offset;
% }
