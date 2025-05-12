# collapse at levels 1 (kingdom) through 7 (species)
for LEVEL in {1..7}; do
  qiime taxa collapse \
    --i-table 250414_StonyCoral_270x200_featureTable.qza \
    --i-taxonomy 250414_270x200_representative-sequences_taxonomy.qza \
    --p-level $LEVEL \
    --o-collapsed-table collapsed-l${LEVEL}.qza
done

